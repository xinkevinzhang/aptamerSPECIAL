import random
import os
from collections import defaultdict

def generate_g_quadruplex(length=40, num_loops=3, loop_length=3):
    """生成G-四链体结构序列"""
    # G-四链体核心结构: G重复序列 + 环序列
    g_repeats = [''.join(['G']*random.randint(3,5)) for _ in range(num_loops+1)]
    loops = [''.join(random.choices('ACU', k=loop_length)) for _ in range(num_loops)]
    
    sequence = g_repeats[0]
    for i in range(num_loops):
        sequence += loops[i] + g_repeats[i+1]
    
    # 如果长度不足，随机添加序列
    if len(sequence) < length:
        sequence += ''.join(random.choices('ACGU', k=length-len(sequence)))
    else:
        sequence = sequence[:length]
    
    return sequence


def generate_pseudoknot(length=40):
    """生成假结结构序列"""
    # 假结结构: 茎环结构 + 交叉配对区域
    stem1_length = random.randint(5,8)
    loop1_length = random.randint(3,5)
    stem2_length = random.randint(4,6)
    loop2_length = random.randint(2,4)
    
    # 生成互补序列
    stem1_part1 = ''.join(random.choices('AUCG', k=stem1_length))
    stem1_part2 = ''.join([{'A':'U','U':'A','C':'G','G':'C'}[c] for c in stem1_part1[::-1]])
    
    stem2_part1 = ''.join(random.choices('AUCG', k=stem2_length))
    stem2_part2 = ''.join([{'A':'U','U':'A','C':'G','G':'C'}[c] for c in stem2_part1[::-1]])
    
    loop1 = ''.join(random.choices('AUCG', k=loop1_length))
    loop2 = ''.join(random.choices('AUCG', k=loop2_length))
    
    # 假结结构排列: stem1_part1 + loop1 + stem2_part1 + loop2 + stem2_part2 + stem1_part2
    sequence = stem1_part1 + loop1 + stem2_part1 + loop2 + stem2_part2 + stem1_part2
    
    # 调整长度
    if len(sequence) < length:
        sequence += ''.join(random.choices('AUCG', k=length-len(sequence)))
    else:
        sequence = sequence[:length]
    
    return sequence


def generate_triplex(length=40):
    """生成三链结构序列"""
    # 三链结构: 双链区域 + 第三链
    duplex_length = random.randint(8,12)
    triplex_length = random.randint(6,10)
    
    # 生成Watson-Crick双链
    strand1 = ''.join(random.choices('ATGC', k=duplex_length))
    strand2 = ''.join([{'A':'T','T':'A','C':'G','G':'C'}[c] for c in strand1])
    
    # 生成Hoogsteen配对的第三链
    triplex = []
    for c in strand1[:triplex_length]:
        if c == 'A':
            triplex.append('T')  # A-T-A三链
        elif c == 'T':
            triplex.append('A')  # T-A-T三链
        elif c == 'G':
            triplex.append('G')  # G-C-G三链
        else:  # C
            triplex.append('C')  # C-G-C三链
    triplex = ''.join(triplex)
    
    # 组合序列
    sequence = strand1 + strand2 + triplex
    
    # 调整长度
    if len(sequence) < length:
        sequence += ''.join(random.choices('ATGC', k=length-len(sequence)))
    else:
        sequence = sequence[:length]
    
    return sequence


def generate_random_sequence(length=40):
    """生成随机序列作为对照"""
    return ''.join(random.choices('ATGC', k=length))


def generate_library(structure_type, count=1000, length=40):
    """生成特定结构类型的序列库"""
    library = []
    label = structure_type.lower().replace(' ', '_')
    
    generators = {
        'g_quadruplex': generate_g_quadruplex,
        'pseudoknot': generate_pseudoknot,
        'triplex': generate_triplex,
        'random': generate_random_sequence
    }
    
    if structure_type not in generators:
        raise ValueError(f"不支持的结构类型: {structure_type}")
    
    for i in range(count):
        seq = generators[structure_type](length=length)
        # 添加标签: {结构类型}_{序号}
        seq_id = f">{label}_seq_{i+1}_{label}"
        library.append(f"{seq_id}\n{seq}")
    
    return library


def main():
    # 创建输出目录
    output_dir = 'aptamer_libraries'
    os.makedirs(output_dir, exist_ok=True)
    
    # 定义要生成的序列类型和数量
    structures = [
        ('g_quadruplex', 2000),
        ('pseudoknot', 2000),
        ('triplex', 2000),
        ('random', 3000)  # 随机序列作为对照
    ]
    
    # 生成并保存各个文库
    all_sequences = []
    for struct_type, count in structures:
        print(f"生成{struct_type}序列: {count}条")
        library = generate_library(struct_type, count=count)
        
        # 保存单独的文库
        filename = os.path.join(output_dir, f'{struct_type}_library.fasta')
        with open(filename, 'w') as f:
            f.write('\n'.join(library))
        
        # 添加到总库
        all_sequences.extend(library)
    
    # 保存合并的文库
    combined_filename = os.path.join(output_dir, 'all_aptamers_combined.fasta')
    with open(combined_filename, 'w') as f:
        f.write('\n'.join(all_sequences))
    
    print(f"序列生成完成!\n单独文库保存在: {output_dir}\n合并文库: {combined_filename}")


if __name__ == '__main__':
    main()