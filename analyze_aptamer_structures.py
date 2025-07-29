import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.cluster import DBSCAN
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
from Levenshtein import distance as edit_distance
from umap import UMAP

# Set matplotlib font for proper English display
plt.rcParams["font.family"] = ["Arial", "sans-serif"]
plt.rcParams["axes.labelsize"] = 12
plt.rcParams["axes.titlesize"] = 14
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10

def parse_structure_file(structure_file):
    """Parse secondary structure file to extract sequence, structure and label information"""
    data = []
    with open(structure_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
        
    for i in range(0, len(lines), 3):
        if i+2 >= len(lines):
            break
        
        seq_id = lines[i]
        sequence = lines[i+1]
        struct_line = lines[i+2]
        structure = struct_line.split()[0]
        mfe = struct_line.split()[1] if len(struct_line.split()) > 1 else 'N/A'
        
        # Extract label from ID
        label = seq_id.split('_')[-1] if '_' in seq_id else 'unknown'
        
        data.append({
            'id': seq_id,
            'sequence': sequence,
            'structure': structure,
            'mfe': mfe,
            'label': label
        })
    
    return pd.DataFrame(data)


def calculate_structure_features(structure):
    """Calculate structural features from structure string"""
    if not structure:
        return {}
    
    # Calculate base pairing statistics
    base_pairs = structure.count('(') + structure.count(')')
    hairpins = structure.count('()')
    
    # Calculate structure complexity
    complexity = len(set(structure)) / len(structure) if structure else 0
    
    # Calculate pairing ratio
    paired_ratio = base_pairs / len(structure) if structure else 0
    
    return {
        'length': len(structure),
        'base_pairs': base_pairs,
        'hairpins': hairpins,
        'complexity': complexity,
        'paired_ratio': paired_ratio,
        'gc_content': (structure.count('G') + structure.count('C')) / len(structure) if structure else 0
    }


def compute_edit_distance_matrix(sequences):
    """Compute edit distance matrix between sequences"""
    n = len(sequences)
    dist_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            dist = edit_distance(sequences[i], sequences[j])
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
    
    return dist_matrix


def cluster_sequences(sequences, eps=5, min_samples=5):
    """Cluster sequences using DBSCAN based on edit distance"""
    # Compute edit distance matrix
    dist_matrix = compute_edit_distance_matrix(sequences)
    
    # Apply DBSCAN clustering
    dbscan = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed')
    clusters = dbscan.fit_predict(dist_matrix)
    
    return clusters


def train_classifier(data, feature_cols, target_col='label'):
    """Train random forest classifier to distinguish structure types"""
    # Prepare features and target
    X = data[feature_cols]
    y = data[target_col]
    
    # Split into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )
    
    # Standardize features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train random forest classifier
    clf = RandomForestClassifier(n_estimators=100, random_state=42)
    clf.fit(X_train_scaled, y_train)
    
    # Evaluate model
    y_pred = clf.predict(X_test_scaled)
    print("Classification Report:")
    print(classification_report(y_test, y_pred))
    
    # Plot confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    plt.figure(figsize=(10, 8))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=clf.classes_, yticklabels=clf.classes_)
    plt.xlabel('Predicted Label')
    plt.ylabel('True Label')
    plt.title('Classification Confusion Matrix')
    plt.savefig('classification_confusion_matrix.png')
    plt.close()
    
    # Get feature importance
    feature_importance = pd.DataFrame({
        'feature': feature_cols,
        'importance': clf.feature_importances_
    }).sort_values('importance', ascending=False)
    
    plt.figure(figsize=(10, 6))
    sns.barplot(x='importance', y='feature', data=feature_importance)
    plt.title('Feature Importance')
    plt.savefig('feature_importance.png')
    plt.close()
    
    return clf, scaler, feature_importance


def visualize_clusters(features, labels, title='Sequence Clusters'):
    """Visualize clustering results using t-SNE"""
    # 保留t-SNE可视化
    tsne = TSNE(n_components=2, random_state=42)
    tsne_results = tsne.fit_transform(features)
    
    plt.figure(figsize=(12, 10))
    scatter = plt.scatter(tsne_results[:, 0], tsne_results[:, 1], c=pd.factorize(labels)[0], cmap='viridis', alpha=0.6)
    plt.legend(handles=scatter.legend_elements()[0], labels=list(pd.unique(labels)))
    plt.title(title)
    # Add axis labels
    plt.xlabel('t-SNE Dimension 1')
    plt.ylabel('t-SNE Dimension 2')
    # Add grid lines
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(f'{title.lower().replace(" ", "_")}.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 添加UMAP可视化
    umap = UMAP(n_components=2, random_state=42)
    umap_results = umap.fit_transform(features)
    
    plt.figure(figsize=(12, 10))
    scatter = plt.scatter(umap_results[:, 0], umap_results[:, 1], c=pd.factorize(labels)[0], cmap='viridis', alpha=0.6)
    plt.legend(handles=scatter.legend_elements()[0], labels=list(pd.unique(labels)))
    plt.title(title)
    # Add axis labels
    plt.xlabel('t-SNE Dimension 1')
    plt.ylabel('t-SNE Dimension 2')
    # Add grid lines
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(f'{title.lower().replace(" ", "_")}_tsne.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 添加UMAP可视化
    umap = UMAP(n_components=2, random_state=42)
    umap_results = umap.fit_transform(features)
    
    plt.figure(figsize=(12, 10))
    scatter = plt.scatter(umap_results[:, 0], umap_results[:, 1], c=pd.factorize(labels)[0], cmap='viridis', alpha=0.6)
    plt.legend(handles=scatter.legend_elements()[0], labels=list(pd.unique(labels)))
    plt.title(title)
    # Add axis labels
    plt.xlabel('t-SNE Dimension 1')
    plt.ylabel('t-SNE Dimension 2')
    # Add grid lines
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(f'{title.lower().replace(" ", "_")}_umap.png', dpi=300, bbox_inches='tight')
    plt.close()


def main():
    # 统一特征列定义
    feature_cols = ['length', 'base_pairs', 'hairpins', 'complexity', 'paired_ratio', 'gc_content']
    
    # 解析结构文件
    print("Loading structural data...")
    df = parse_structure_file('secondary_structures.txt')
    
    # Calculate structural features
    print("Calculating structural features...")
    struct_features = df['structure'].apply(calculate_structure_features)
    struct_df = pd.DataFrame(struct_features.tolist())
    df = pd.concat([df, struct_df], axis=1)
    
    # Save feature data
    df.to_csv('aptamer_features.csv', index=False)
    print("Feature data saved to aptamer_features.csv")
    
    # 新增：生成分析报告
    with open('analysis_report.txt', 'w') as f:
        f.write("Aptamer Structure Analysis Report\n")
        f.write("================================\n\n")
        f.write(f"Total sequences analyzed: {len(df)}\n")
        f.write(f"Structure types distribution:\n")
        f.write(f"{df['label'].value_counts().to_string()}\n\n")
        f.write(f"Feature statistics:\n")
        f.write(f"{df[feature_cols].describe().to_string()}\n")
    print("Analysis report generated: analysis_report.txt")
    
    # 样本选择代码 - 此行是正确的分层抽样
    sample_df = df.groupby('label', group_keys=False).apply(lambda x: x.sample(min(100, len(x))))
    
    # 聚类分析
    print("Performing sequence clustering...")
    # 使用已有的分层抽样结果，删除重复的随机抽样代码
    clusters = cluster_sequences(sample_df['sequence'].tolist())
    sample_df['cluster'] = clusters
    
    # 可视化聚类结果前添加
    print(f"Visualizing {len(sample_df)} samples for cluster analysis")
    print(f"Features used: {feature_cols}")
    visualize_clusters(
        sample_df[feature_cols].values,
        sample_df['label'],
        title='Sequence Type Clusters'
    )
    
    visualize_clusters(
        sample_df[['length', 'base_pairs', 'hairpins', 'complexity', 'paired_ratio']].values,
        sample_df['cluster'],
        title='Edit Distance Clusters'
    )
    
    # New: Feature distribution box plots
    print("Generating feature distribution plots...")
    feature_cols = ['length', 'base_pairs', 'hairpins', 'complexity', 'paired_ratio', 'gc_content']
    plt.figure(figsize=(15, 10))
    for i, col in enumerate(feature_cols):
        plt.subplot(2, 3, i+1)
        sns.boxplot(x='label', y=col, data=df)
        plt.title(f'Distribution of {col} by Structure Type')
        plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('feature_distributions_by_type.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # New: Feature correlation heatmap
    print("Generating feature correlation heatmap...")
    corr = df[feature_cols].corr()
    plt.figure(figsize=(10, 8))
    sns.heatmap(corr, annot=True, cmap='coolwarm', vmin=-1, vmax=1)
    plt.title('Feature Correlation Heatmap')
    plt.savefig('feature_correlation_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Train classifier
    print("Training structure type classifier...")
    # 修复前：多处定义feature_cols
    feature_cols = ['length', 'base_pairs', 'hairpins', 'complexity', 'paired_ratio', 'gc_content']
    
    # Filter labeled data
    labeled_df = df[df['label'] != 'unknown']
    if len(labeled_df) < 10:
        print("Warning: Insufficient labeled data for classifier training")
    else:
        clf, scaler, feature_importance = train_classifier(labeled_df, feature_cols)
        print("Classifier training complete, feature importance:")
        print(feature_importance)
        
        # Predict labels for all samples
        df['predicted_label'] = clf.predict(scaler.transform(df[feature_cols]))
        
        # Save data with predictions
        df.to_csv('aptamer_with_predictions.csv', index=False)
        print("Data with predictions saved to aptamer_with_predictions.csv")

if __name__ == '__main__':
    main()