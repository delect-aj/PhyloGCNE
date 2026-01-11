#!/bin/python

import json
import logging
import shutil
import argparse
import os
import time
import random
import csv
from typing import List
from tqdm import tqdm
from typing import Dict, List
from collections import OrderedDict

import biom
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skbio import TreeNode

from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, roc_auc_score

import torch
import torch.nn as nn
from torch import Tensor
import torch.optim as optim
import torch.backends.cudnn as cudnn
import torch_geometric
import torch.nn.functional as F
from torch_geometric.loader import DataLoader
from torch_geometric.data import Dataset, Data
from torch_geometric.nn.conv import MessagePassing
from torch_geometric.utils import add_self_loops, degree
from torch_geometric.nn import global_mean_pool, global_max_pool
from torch_geometric.typing import Adj
import torch_geometric.graphgym.register as register


def random_seed(seed=42, rank=0):
    """Sets the random seed for reproducibility."""
    torch.manual_seed(seed + rank)
    np.random.seed(seed + rank)
    random.seed(seed + rank)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed + rank)
    logging.info(f"Random seed set to {seed + rank}")


def get_logger(filepath: str, level=logging.INFO, name: str = None) -> logging.Logger:
    """Configures and returns a logger."""
    logger_name = name if name else os.path.splitext(os.path.basename(filepath))[0]
    logger = logging.getLogger(logger_name)
    logger.setLevel(level)
    if not logger.handlers:
        file_handler = logging.FileHandler(filepath, mode='w')
        file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
        
        stream_handler = logging.StreamHandler()
        stream_formatter = logging.Formatter('%(levelname)s: %(message)s')
        stream_handler.setFormatter(stream_formatter)
        logger.addHandler(stream_handler)
    return logger

  
def gaussian_kernel(d, rho=2):
    """
    weight = exp(-2 * ρ * d^2)
    """
    return np.exp(-2 * rho * d ** 2)


class MicroGraphAblationDataset(Dataset):
    """
    """
    def __init__(self,
                 biom_file_path: str, 
                 metadata_file_path: str,
                 phylogeny_file_path: str,
                 label_column: str = 'group', transform=None):
                    
        super().__init__(transform=transform)

        self.table = biom.load_table(biom_file_path)
        full_tree = TreeNode.read(phylogeny_file_path)
        self.fid = self.table.ids(axis='observation')
        df_metadata = pd.read_csv(metadata_file_path, index_col=0, sep='\t')
        # self.table = self.table.norm(axis='sample', inplace=False)
        self.table = self.table.to_dataframe().T
        
        data = {}
        self.all_nodes_in_tree = [node.name for node in full_tree.postorder()]
        
        self.data_list = []
        for sample_id in tqdm(self.table.index.values):
            sample_data = np.array(self.table.loc[sample_id])
            label = df_metadata.loc[sample_id, label_column]
            
            data = self._build_graph_from_sample(self.fid, sample_data, label, sample_id, full_tree)
            self.data_list.append(data)
        

    def _build_graph_from_sample(self, otu_list, abundance_list, label, sample_id, 
                                 full_tree):
        leaf_abundance_map = {otu: abund for otu, abund in zip(otu_list, abundance_list)}
        
        tree_sample = full_tree
        all_nodes_in_tree = list(tree_sample.traverse(self_before=True))
        node_to_idx = {node.name: i for i, node in enumerate(all_nodes_in_tree)}

        full_abundance_map = {}
        x = torch.zeros((len(all_nodes_in_tree), 1), dtype=torch.float)
        edge_list, edge_weight_list = [], []
        node_id = []
        for node in tree_sample.postorder():
            node_id.append(node_to_idx[node.name])
            if node.is_tip():
                full_abundance_map[node.name] = leaf_abundance_map.get(node.name, 0.0)
                x[node_to_idx[node.name]] = full_abundance_map[node.name]
            else:
                child_abundances = [full_abundance_map.get(child.name, 0.0) for child in node.children]
                full_abundance_map[node.name] = sum(child_abundances)
                x[node_to_idx[node.name]] = full_abundance_map[node.name]
                
            if not node.is_root():
                node_idx = node_to_idx[node.name]
                parent_idx = node_to_idx[node.parent.name]
                parent_idx = node_to_idx[node.parent.name]
                edge_list.extend([[node_idx, parent_idx]])
                distance = node.length if node.length is not None else 0.0
                weight = gaussian_kernel(distance, rho=2) 
                edge_weight_list.extend([weight])

        edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous()
        y = torch.tensor([label], dtype=torch.long)
        node_id = torch.tensor(node_id, dtype=torch.long)
        edge_attr = torch.tensor(edge_weight_list, dtype=torch.float).view(-1, 1)
        
        data = Data(x=x, node_id=node_id, edge_index=edge_index, edge_attr=edge_attr, y=y, sample_id=sample_id)

        return data

    def len(self) -> int:
        return self.table.shape[0]

    def get(self, idx: int):
            
        return self.data_list[idx]
    
    def collate_fn(self, batch):
        return torch_geometric.data.Batch.from_data_list(batch)
      
    def __call__(self):
        """Get num nodes"""
        return len(self.all_nodes_in_tree)


class _FFNBlock(nn.Module):
    def __init__(self, dim, dim_ff, dropout, act_fn, norm_type):
        super().__init__()
        self.linear1 = nn.Linear(dim, dim_ff)
        self.act_fn = act_fn
        self.dropout1 = nn.Dropout(dropout)
        self.linear2 = nn.Linear(dim_ff, dim)
        self.dropout2 = nn.Dropout(dropout)
        
        self.norm = None
        if norm_type == 'batch_norm':
            self.norm = nn.BatchNorm1d(dim)
        elif norm_type == 'layer_norm':
            self.norm = nn.LayerNorm(dim)
            
    def forward(self, x):
        identity = x
        out = self.dropout1(self.act_fn(self.linear1(x)))
        out = self.dropout2(self.linear2(out))
        out = identity + out
        if self.norm:
            out = self.norm(out)
        return out
        

class GCNConvWithEdges(MessagePassing):
    def __init__(self, in_channels: int, out_channels: int, edge_dim: int, 
                 add_self_loops: bool = True, bias: bool = True):
        super().__init__(aggr='add')

        self.add_self_loops_flag = add_self_loops
        self.in_channels = in_channels
        self.out_channels = out_channels

        self.lin = nn.Linear(in_channels, out_channels, bias=False)
        self.edge_lin = nn.Linear(edge_dim, out_channels, bias=False)

        if bias:
            self.bias = nn.Parameter(torch.Tensor(out_channels))
        else:
            self.register_parameter('bias', None)

        self.reset_parameters()

    def reset_parameters(self):
        nn.init.xavier_uniform_(self.lin.weight)
        nn.init.xavier_uniform_(self.edge_lin.weight)
        if self.bias is not None:
            nn.init.zeros_(self.bias)

    def forward(self, x: Tensor, edge_index: Tensor, edge_attr: Tensor) -> Tensor:
        num_nodes = x.size(0)

        if self.add_self_loops_flag:
            edge_index, edge_attr = self.add_self_loops_with_attr(
                edge_index, edge_attr, num_nodes
            )

        row, col = edge_index
        deg = degree(col, num_nodes, dtype=x.dtype)
        deg_inv_sqrt = deg.pow(-0.5)
        deg_inv_sqrt[deg_inv_sqrt == float('inf')] = 0
        norm = deg_inv_sqrt[row] * deg_inv_sqrt[col]

        x = self.lin(x)                  
        edge_embedding = self.edge_lin(edge_attr) 

        out = self.propagate(edge_index, x=x, edge_attr=edge_embedding, norm=norm)

        if self.bias is not None:
            out += self.bias

        return out

    def message(self, x_j: Tensor, edge_attr: Tensor, norm: Tensor) -> Tensor:
        return norm.view(-1, 1) * (x_j + edge_attr)

    def add_self_loops_with_attr(self, edge_index, edge_attr, num_nodes):
        edge_index, _ = add_self_loops(edge_index, num_nodes=num_nodes)
        loop_attr = torch.zeros((num_nodes, edge_attr.size(1)), 
                                device=edge_attr.device, dtype=edge_attr.dtype)
        edge_attr = torch.cat([edge_attr, loop_attr], dim=0)
        return edge_index, edge_attr


class BaseEdgeGNNLayer(nn.Module):
    def __init__(self, dim_in, dim_out, dropout, residual, ffn, norm_type, act):
        super().__init__()
        self.dim_in = dim_in
        self.dim_out = dim_out
        self.dropout = dropout
        self.residual = residual
        self.norm_type = norm_type
        self.ffn = ffn
        
        self.conv = None

        self.norm = None
        if self.norm_type == 'batch_norm':
            self.norm = nn.BatchNorm1d(dim_out)
        elif self.norm_type == 'layer_norm':
            self.norm = nn.LayerNorm(dim_out)
        
        self.act = nn.Sequential(
            register.act_dict[act](),
            nn.Dropout(self.dropout),
        )

        if self.residual and self.dim_in != self.dim_out:
            self.residual_lin = nn.Linear(dim_in, dim_out)
        else:
            self.residual_lin = nn.Identity()

        if self.ffn:
            self.ffn_block = _FFNBlock(
                dim=dim_out, 
                dim_ff=dim_out * 2, 
                dropout=dropout, 
                act_fn=register.act_dict[act](), 
                norm_type=norm_type
            )

    def forward(self, batch):
        x_in = batch.x
        x_out = self.conv(batch.x, batch.edge_index, batch.edge_attr)
        if self.norm:
            x_out = self.norm(x_out)
        x_out = self.act(x_out)
        if self.residual:
            x_out = self.residual_lin(x_in) + x_out
        if self.ffn:
            x_out = self.ffn_block(x_out)
        batch.x = x_out
        return batch


class GCNConvEdgeLayer(BaseEdgeGNNLayer):
    def __init__(self, dim_in, dim_out, edge_dim, dropout, residual, ffn,
                 norm_type, act='relu'):
        super().__init__(dim_in, dim_out, dropout, residual, ffn, norm_type, act)
        self.conv = GCNConvWithEdges(dim_in, dim_out, edge_dim=edge_dim, bias=True)


class GCNEClassifier(nn.Module):
    """
    """
    def __init__(self, 
                 num_nodes: int,
                 in_channels: int, 
                 hidden_channels: int, 
                 num_classes: int = 2,
                 edge_dim: int = 1,
                 num_layers: int = 3,
                 dropout_rate: float = 0.5, 
                 residual: bool = True,
                 ffn: bool = False,
                 norm_type: str = None,
                 act: str = 'relu',
                 readout: str = 'mean'):
        super().__init__()
        
        self.embedding = nn.Embedding(num_nodes, in_channels)
        self.readout = readout
        self.num_classes = num_classes
        self.conv_layers = nn.ModuleList()
        self.linear = nn.Linear(1, in_channels)
        
        self.conv_layers.append(
            GCNConvEdgeLayer(in_channels, hidden_channels, edge_dim, dropout_rate, 
                             residual=residual, ffn=ffn, norm_type=norm_type, act=act)
        )
        
        for _ in range(num_layers - 1):
            self.conv_layers.append(
                GCNConvEdgeLayer(hidden_channels, hidden_channels, edge_dim, dropout_rate, 
                                 residual=residual, ffn=ffn, norm_type=norm_type, act=act)
            )
        
        if self.readout == 'attention':
            self.attention = AttentionReadout(hidden_channels)
        

        output_dim = num_classes if num_classes > 2 else 1
        pool_dim = hidden_channels * 2 if self.readout == 'mean_max' else hidden_channels

        self.head = nn.Sequential(
            nn.Linear(pool_dim, hidden_channels // 2),
            register.act_dict[act](),  
            nn.Dropout(p=dropout_rate),
            nn.Linear(hidden_channels // 2, output_dim)
        )
        
    def forward(self, data):
        data.x = self.linear(data.x) + self.embedding(data.node_id)
        
        for layer in self.conv_layers:
            data = layer(data)
            
        if self.readout == 'mean':
            x = global_mean_pool(data.x, data.batch)
        elif self.readout == 'max':
            x = global_max_pool(data.x, data.batch)
        elif self.readout == 'attention':
            x, attn_score = self.attention(data.x, data.batch)
        elif self.readout == 'mean_max':
            x = torch.cat([global_mean_pool(data.x, data.batch), global_max_pool(data.x, data.batch)], dim=1)
        
        x = self.head(x)
        
        if self.num_classes == 2:
            return x.squeeze(-1)
        
        return x
  

def train_epoch(model, loader, optimizer, device, logger, epoch_num, scheduler, args):
    model.train()
    losses_m = AverageMeter('Loss', ':.4e')
    
    all_probs = []
    all_preds = []
    all_labels = []
    all_train_sample_ids = []

    for batch_idx, batch_data in enumerate(loader):
        batch_data = batch_data.to(device, non_blocking=True)
        y = batch_data.y 

        if hasattr(batch_data, 'sample_id'):
            if isinstance(batch_data.sample_id, list):
                all_train_sample_ids.extend(batch_data.sample_id)
            else:
                all_train_sample_ids.append(batch_data.sample_id)
        else:
            if batch_idx == 0:
                 logger.warning("Train loop: 'sample_id' attribute not found.")

        optimizer.zero_grad()
        logits = model(batch_data)
        
        is_multiclass = logits.dim() > 1 and logits.size(1) > 1
        
        if is_multiclass:
            criterion = nn.CrossEntropyLoss()
            loss = criterion(logits, y.long())
            
            probs = torch.softmax(logits, dim=1) 
            preds = torch.argmax(probs, dim=1)   
            
        else:
            criterion = nn.BCEWithLogitsLoss()
            loss = criterion(logits.view(-1), y.float().view(-1))
            
            probs = torch.sigmoid(logits).view(-1)
            preds = (probs > 0.5).long()          

        loss.backward()
        
        if hasattr(args, 'max_grad_norm') and args.max_grad_norm > 0:
            torch.nn.utils.clip_grad_norm_(model.parameters(), args.max_grad_norm)
            
        optimizer.step()
        if scheduler:
            scheduler.step()

        losses_m.update(loss.item(), batch_data.num_graphs)
        
        all_probs.append(probs.detach().cpu().numpy())
        all_preds.append(preds.detach().cpu().numpy())
        all_labels.append(y.detach().cpu().numpy())

    all_labels = np.concatenate(all_labels)
    all_preds = np.concatenate(all_preds)
    all_probs = np.concatenate(all_probs)

    acc = accuracy_score(all_labels, all_preds)
    
    is_multiclass_metric = (all_probs.ndim == 2) and (all_probs.shape[1] > 1)

    if is_multiclass_metric:
        f1 = f1_score(all_labels, all_preds, average='macro', zero_division=0)
        try:
            auc_score = roc_auc_score(all_labels, all_probs, multi_class='ovr', average='macro')
        except ValueError:
            auc_score = 0.0
    else:
        f1 = f1_score(all_labels, all_preds, zero_division=0)
        try:
            auc_score = roc_auc_score(all_labels, all_probs)
        except ValueError:
            auc_score = 0.0

    metrics = OrderedDict([('loss', losses_m.avg), ('acc', acc), ('f1', f1), ('auc', auc_score)])
    
    return metrics, all_train_sample_ids


def evaluate_epoch(model, loader, device, logger, epoch_num_str="Eval"):
    model.eval()
    losses_m = AverageMeter('Loss', ':.4e')
    
    all_probs = []
    all_preds = []
    all_labels = []
    all_sample_ids = []

    with torch.no_grad():
        for batch_data in loader:
            batch_data = batch_data.to(device, non_blocking=True)
            y = batch_data.y 
            
            logits = model(batch_data)
            is_multiclass = logits.dim() > 1 and logits.size(1) > 1

            if is_multiclass:
                criterion = nn.CrossEntropyLoss()
                loss = criterion(logits, y.long()) 
                
                probs = torch.softmax(logits, dim=1) # (N, C)
                preds = torch.argmax(probs, dim=1)   # (N)
            else:
                criterion = nn.BCEWithLogitsLoss()
                loss = criterion(logits.view(-1), y.float().view(-1)) 
                
                probs = torch.sigmoid(logits).view(-1)
                preds = (probs > 0.5).long()

            if loss is not None:
                losses_m.update(loss.item(), batch_data.num_graphs)

            all_probs.append(probs.cpu().numpy())
            all_preds.append(preds.cpu().numpy())
            all_labels.append(y.cpu().numpy())

            if hasattr(batch_data, 'sample_id'):
                if isinstance(batch_data.sample_id, list):
                    all_sample_ids.extend(batch_data.sample_id)
                else:
                    all_sample_ids.append(batch_data.sample_id)

    all_labels = np.concatenate(all_labels)
    all_preds = np.concatenate(all_preds)
    all_probs = np.concatenate(all_probs)
    
    acc = accuracy_score(all_labels, all_preds)

    is_multiclass_metric = (all_probs.ndim == 2) and (all_probs.shape[1] > 1)

    if is_multiclass_metric:
        avg_method = 'macro' 
        f1 = f1_score(all_labels, all_preds, average=avg_method, zero_division=0)
        prec = precision_score(all_labels, all_preds, average=avg_method, zero_division=0)
        rec = recall_score(all_labels, all_preds, average=avg_method, zero_division=0)
        try:
            auc_score = roc_auc_score(all_labels, all_probs, multi_class='ovr', average=avg_method)
        except ValueError:
            auc_score = 0.0
    else:
        f1 = f1_score(all_labels, all_preds, zero_division=0)
        prec = precision_score(all_labels, all_preds, zero_division=0)
        rec = recall_score(all_labels, all_preds, zero_division=0)
        try:
            auc_score = roc_auc_score(all_labels, all_probs)
        except ValueError:
            auc_score = 0.0

    metrics = OrderedDict([
        ('loss', losses_m.avg), 
        ('acc', acc), 
        ('f1', f1), 
        ('auc', auc_score), 
        ('precision', prec), 
        ('recall', rec)
    ])
    
    if len(all_sample_ids) != len(all_labels):
        logger.warning(f"Sample ID count ({len(all_sample_ids)}) does not match label count ({len(all_labels)}).")
        all_sample_ids = None
    
    return metrics, all_labels, all_probs, all_sample_ids



def calculate_node_importance_beta(model, loader, device, all_nodes_names, num_classes, output_dir, fold, 
                                   n_samples=50, stdev_spread=0.15, beta_alpha=2.0, beta_beta=2.0):
    """
    Computing Feature Importance Using Beta Distribution Noise (BAS-Grad Style)
    Args:
        beta_alpha, beta_beta: Shape parameters of the Beta distribution.
        
        (2, 2) resembles a bell-shaped distribution (bounded, peaked in the middle).
        
        (1, 1) corresponds to a uniform distribution.
        
        (0.5, 0.5) yields a U-shaped distribution.
    """
    model.eval()
    num_tree_nodes = len(all_nodes_names)
    
    node_grad_sum = torch.zeros(num_tree_nodes, device=device)
    node_counts = torch.zeros(num_tree_nodes, device=device)

    print(f"Computing Feature Importance with Beta Noise (alpha={beta_alpha}, beta={beta_beta})...")
    
    for param in model.parameters():
        param.requires_grad = False

    beta_dist = torch.distributions.Beta(
        torch.tensor([beta_alpha], device=device), 
        torch.tensor([beta_beta], device=device)
    )

    for batch_idx, batch_data in enumerate(tqdm(loader, desc="Calculating")):
        batch_data = batch_data.to(device)
        original_x = batch_data.x.clone().detach()
        
        std_per_feature = original_x.std(dim=0, keepdim=True) 
        std_per_feature[std_per_feature == 0] = 1e-5
        sigma = std_per_feature * stdev_spread

        batch_grads_sum = torch.zeros_like(original_x)

        for i in range(n_samples):
            beta_noise = beta_dist.sample(original_x.shape).squeeze(-1)
            centered_noise = (beta_noise - 0.5)

            noise = centered_noise * sigma 

            curr_noisy_x = torch.clamp(original_x + noise, min=0.0)
            
            curr_noisy_x.requires_grad = True
            batch_data.x = curr_noisy_x
            
            logits = model(batch_data)
            
            if num_classes > 2:
                score, _ = torch.max(logits, dim=1)
            else:
                score = logits.squeeze()
            
            model.zero_grad()
            score.sum().backward()
            
            if curr_noisy_x.grad is not None:
                batch_grads_sum += curr_noisy_x.grad.detach()
            
            curr_noisy_x.grad = None
            
        avg_grads = batch_grads_sum / n_samples
        
        if avg_grads.dim() > 1:
            avg_grads = avg_grads.sum(dim=1)

        avg_grads = avg_grads.view(-1)
        node_ids = batch_data.node_id
        
        node_grad_sum.scatter_add_(0, node_ids, avg_grads)
        node_counts.scatter_add_(0, node_ids, torch.ones_like(avg_grads))
        
        batch_data.x = original_x

    for param in model.parameters():
        param.requires_grad = True

    node_counts[node_counts == 0] = 1.0 
    
    avg_signed_saliency = node_grad_sum / node_counts
    
    importance_scores = avg_signed_saliency.abs().cpu().numpy()
    directions = avg_signed_saliency.cpu().numpy()
    
    df = pd.DataFrame({
        'Node_Name': all_nodes_names,
        'Importance_Score': importance_scores,
        'Direction': directions})
    
    df = df.sort_values(by='Importance_Score', ascending=False)
    
    filename = f"feature_importance_{fold}.csv"
    save_path = os.path.join(output_dir, filename)
    df.to_csv(save_path, index=False)
    print(f"Feature importance (Beta Noise) successfully saved to: {save_path}")
    
      
def plot_metrics(history: Dict[str, List[float]], output_dir: str, fold_num: str, logger: logging.Logger):
    """Plots training and test metrics (Loss and AUC) over epochs and saves the figure."""
    if not history or not history['train_loss']:
        logger.warning("History is empty, skipping plotting.")
        return

    num_epochs = len(history['train_loss'])
    epochs_range = range(1, num_epochs + 1)

    plt.style.use('seaborn-v0_8-whitegrid')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3))
    # Plot 1: Loss vs. Epochs
    ax1.plot(epochs_range, history['train_loss'], 'o-', label='Train Loss', color='blue', markersize=2)
    ax1.plot(epochs_range, history['test_loss'], 'o-', label='Test Loss', color='red', markersize=2)
    ax1.set_title('Training and Test Loss', fontsize=10)
    ax1.set_xlabel('Epoch', fontsize=12)
    ax1.set_ylabel('Loss', fontsize=12)
    ax1.legend()
    ax1.tick_params(axis='both', which='major', labelsize=6)

    # Plot 2: AUC vs. Epochs
    ax2.plot(epochs_range, history['train_auc'], 'o-', label='Train AUC', color='blue', markersize=2)
    ax2.plot(epochs_range, history['test_auc'], 'o-', label='Test AUC', color='red', markersize=2)
    ax2.set_title('Training and Test AUC', fontsize=10)
    ax2.set_xlabel('Epoch', fontsize=12)
    ax2.set_ylabel('AUC', fontsize=12)
    ax2.legend()
    ax2.tick_params(axis='both', which='major', labelsize=6)

    fig.suptitle(f'Metrics for Fold: {fold_num}', fontsize=12, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    plot_path = os.path.join(output_dir, f'metrics_plot_{fold_num}.png')
    try:
        plt.savefig(plot_path, dpi=300)
        logger.info(f"Metrics plot saved to {plot_path}")
    except Exception as e:
        logger.error(f"Failed to save metrics plot to '{plot_path}': {e}", exc_info=True)
    plt.close()


class AverageMeter(object):
    """Computes and stores the average and current value."""
    def __init__(self, name, fmt=':f'):
        self.name, self.fmt = name, fmt
        self.reset()
    def reset(self):
        self.val, self.avg, self.sum, self.count = 0, 0, 0, 0
    def update(self, val, n=1):
        self.val, self.sum, self.count = val, self.sum + val * n, self.count + n
        self.avg = self.sum / self.count if self.count > 0 else 0
        

def main():
    parser = argparse.ArgumentParser(description="Microbe GCN Model Training")
    
    # Path Arguments
    parser.add_argument("--train_table", type=str,  help="Data directory containing train data.")
    parser.add_argument("--test_table", type=str, help="Data directory containing test data.")
    parser.add_argument("--metadata_filename", type=str, help="Metadata filename.")
    parser.add_argument("--phylogeny_file_path", type=str, help="Path to the phylogeny tree file.")
    parser.add_argument("--fold", type=str, help="Current test study/fold name for LODO evaluation.")
    parser.add_argument("--output_dir", type=str, default="./graph_model_crc_workdir/v19", help="Root output directory.")
    parser.add_argument("--feature_importance", type=bool, default=False, help="")

    # Data Arguments
    parser.add_argument("--label_column", type=str, default="group", help="Name of the label column in the metadata file.")

    # Model Parameters
    parser.add_argument("--hidden_channels", type=int, default=100, help="Number of hidden units in GCN layers.")
    parser.add_argument("--heads", type=int, default=1, help="Number of attention heads in the GATConv layer.")
    parser.add_argument("--dropout_rate", type=float, default=0.0, help="Dropout rate used in the classifier head.")
    parser.add_argument("--num_layers", type=int, default=5, help="Total number of GCN layers in the model.")
    parser.add_argument("--act", type=str, default='relu', help="Activation function used in FFN layers.")
    parser.add_argument("--ffn", type=bool, default=True, help="Enable Feed-Forward Network (FFN) blocks in GCN layers.")
    parser.add_argument("--residual", type=bool, default=True, help="Enable residual connections in GCN layers.")
    parser.add_argument("--norm_type", type=str, default='batch_norm', help="Type of normalization to use in GCN layers.")
    parser.add_argument("--readout", type=str, default='mean', help="Readout operation to use in GCN layers.")

    # Training Hyperparameters
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch_size", type=int, default=16)
    parser.add_argument("--learning_rate", type=float, default=5e-4)
    parser.add_argument("--weight_decay", type=float, default=0.0)
    parser.add_argument("--optimizer", type=str, default="adamw", choices=["adamw", "adam"])
    parser.add_argument("--lr_scheduler_type", type=str, default="none", choices=['linear', 'cosine', 'none'])
    parser.add_argument("--lr_warmup_ratio", type=float, default=0.0)
    parser.add_argument("--max_grad_norm", type=float, default=1.0, help="Maximum gradient norm for clipping.")

    # Other Parameters
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--stdev_spread", type=float, default=0.15)
    parser.add_argument("--beta_alpha", type=float, default=2)
    parser.add_argument("--beta_beta", type=float, default=2)
    parser.add_argument("--device", type=str, default="cuda")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--log_interval", type=int, default=20)

    args = parser.parse_args()

    # --- Setup Output Dir ---
    os.makedirs(args.output_dir, exist_ok=True)

    config_flags = []
    if args.ffn: config_flags.append('F')       # F for FFN
    if args.residual: config_flags.append('R')   # R for Residual
    if args.norm_type == 'batch_norm': config_flags.append('B') # BN for BatchNorm
    if args.norm_type == 'layer_norm': config_flags.append('L') # LN for LayerNorm
    config_str = "".join(config_flags)

    # Device and Seed
    logger = get_logger(os.path.join(args.output_dir, f"train_{args.fold}.log"), name=f"TrainFold{args.fold}")
    random_seed(args.seed)
    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    if device.type == 'cuda': cudnn.benchmark = True
    logger.info(f"Using device: {device}")

    # --- Pre-check number of classes ---
    try:
        df_meta = pd.read_csv(args.metadata_filename, sep='\t')
        unique_labels = df_meta[args.label_column].unique()
        num_classes = len(unique_labels)
        logger.info(f"Detected {num_classes} classes in metadata column '{args.label_column}': {unique_labels}")
    except Exception as e:
        logger.warning(f"Could not determine num_classes from metadata file: {e}. Defaulting to 2.")
        num_classes = 2

    # Create Datasets and DataLoaders
    logger.info("Loading datasets...")

    train_dataset = MicroGraphAblationDataset(biom_file_path=args.train_table, metadata_file_path=args.metadata_filename,
        phylogeny_file_path=args.phylogeny_file_path, label_column=args.label_column
    )

    test_dataset = MicroGraphAblationDataset(biom_file_path=args.test_table, metadata_file_path=args.metadata_filename,
        phylogeny_file_path=args.phylogeny_file_path, label_column=args.label_column
    )
    
    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True, num_workers=args.workers)
    test_loader = DataLoader(test_dataset, batch_size=args.batch_size, shuffle=False, num_workers=args.workers)
    logger.info(f"Datasets loaded. Train samples: {len(train_dataset)}, Test samples: {len(test_dataset)}")

    # Initialize model
    num_nodes = train_dataset()
    logger.info(f"Number of nodes: {num_nodes}")
    in_channels = 32
    logger.info(f"Detected input node feature dimension: {in_channels}")
    
    model = GCNEClassifier(
                in_channels=in_channels,
                hidden_channels=args.hidden_channels,
                num_nodes=num_nodes,
                num_classes=num_classes, 
                num_layers=args.num_layers,
                dropout_rate=args.dropout_rate,
                residual=args.residual,
                ffn=args.ffn,
                norm_type=args.norm_type,
                act=args.act,
                readout=args.readout
            ).to(device)
    logger.info(f"Model initialized. Trainable parameters: {sum(p.numel() for p in model.parameters() if p.requires_grad):,}")
    best_model_path = os.path.join(args.output_dir, f"{args.fold}_test_model.pt")
    if args.feature_importance == False:
        # Optimizer, and Scheduler
        optimizer = optim.AdamW(model.parameters(), lr=args.learning_rate, weight_decay=args.weight_decay)
        
        scheduler = None
        if args.lr_scheduler_type != 'none':
            try:
                from transformers import get_scheduler
                total_steps = len(train_loader) * args.epochs
                warmup_steps = int(total_steps * args.lr_warmup_ratio)
                scheduler = get_scheduler(args.lr_scheduler_type, optimizer, num_warmup_steps=warmup_steps, num_training_steps=total_steps)
                logger.info(f"Using '{args.lr_scheduler_type}' scheduler with {warmup_steps} warmup steps.")
            except ImportError:
                logger.warning("`transformers` library not found. Skipping advanced LR scheduler.")
        
        # Training Loop
        logger.info(f"Starting training for fold {args.fold}...")
        best_test_auc = -1.0
        epochs_no_improve = 0
        
        history = {
            'train_loss': [], 'test_loss': [],
            'train_auc': [], 'test_auc': []
        }
        
        for epoch in range(1, args.epochs + 1):
            train_metrics, _ = train_epoch(model, train_loader, optimizer, device, logger, epoch, scheduler, args)
            test_metrics, _, _, _ = evaluate_epoch(model, test_loader, device, logger, f"Test Ep {epoch}")
    
            history['train_loss'].append(train_metrics['loss'])
            history['train_auc'].append(train_metrics['auc'])
            history['test_loss'].append(test_metrics['loss'])
            history['test_auc'].append(test_metrics['auc'])
    
            # Save best model
            is_best = test_metrics['auc'] > best_test_auc
            if is_best:
                best_test_auc = test_metrics['auc']
                epochs_no_improve = 0
                torch.save(model.state_dict(), best_model_path)
            else:
                epochs_no_improve += 1
                
        logger.info(f"Training for fold {args.fold} finished. Best Test AUC: {best_test_auc:.4f}")
        plot_metrics(history, args.output_dir, args.fold, logger)
    
        # --- Final: Load Best Model and Save Predictions ---
        logger.info(f"Loading best model from {best_model_path} for final prediction...")
        model.load_state_dict(torch.load(best_model_path))
        model.eval()
    
        # Re-run evaluation on test set using the best model
        metrics_final, all_labels_final, all_probs_final, all_sample_ids_final = evaluate_epoch(
            model, test_loader, device, logger, epoch_num_str="FinalEval"
        )
    
        # Prepare DataFrame for saving
        prediction_file_path = os.path.join(args.output_dir, f"predictions_{args.fold}.csv")
        
        if all_sample_ids_final is None:
            logger.warning("Sample IDs were lost during evaluation. Creating dummy IDs.")
            all_sample_ids_final = [f"sample_{i}" for i in range(len(all_labels_final))]
    
        if num_classes > 2:
            data_dict = {'SampleID': all_sample_ids_final}
            for i in range(num_classes):
                data_dict[f'Prob_Class_{i}'] = all_probs_final[:, i]
            data_dict['True_Label'] = all_labels_final
            
            df_preds = pd.DataFrame(data_dict)
            logger.info(f"Saving multi-class predictions to {prediction_file_path}")
            
        else:
            df_preds = pd.DataFrame({
                'SampleID': all_sample_ids_final,
                'Prob_Class_1': all_probs_final,
                'True_Label': all_labels_final
            })
            logger.info(f"Saving binary classification predictions to {prediction_file_path}")
    
        df_preds.to_csv(prediction_file_path, index=False)
        logger.info("Predictions saved successfully.")
    else: 
        ### model explain
        all_node_names = train_dataset.all_nodes_in_tree
        model.load_state_dict(torch.load(best_model_path))
        calculate_node_importance_beta(
                model=model, 
                loader=test_loader,  
                device=device, 
                all_nodes_names=all_node_names, 
                num_classes=num_classes,
                output_dir=args.output_dir,
                fold=args.fold,
                stdev_spread=args.stdev_spread,
                beta_alpha=args.beta_alpha,
                beta_beta=args.beta_beta
            )
    
if __name__ == "__main__":
    main()
