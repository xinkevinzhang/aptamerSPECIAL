import argparse
import subprocess
import os

def predict_secondary_structures(input_file, output_file, rnafold_path):
    # 确保输出目录存在
    dir_name = os.path.dirname(output_file)
    if dir_name:
        os.makedirs(dir_name, exist_ok=True)
    
    # 直接构建RNAfold命令
    command = [rnafold_path, '--noPS', input_file]
    
    # 执行命令并捕获输出
    result = subprocess.run(
        command,
        capture_output=True,
        text=True
    )
    
    # 保存结果
    with open(output_file, 'w') as f:
        f.write(result.stdout)
    
    print(f"预测完成，结果已保存至: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='预测RNA二级结构')
    parser.add_argument('input', help='输入FASTA文件路径')
    parser.add_argument('-o', '--output', default='secondary_structures.txt', help='输出文件路径')
    parser.add_argument('--rnafold', 
                        default='D:\Program Files (x86)\ViennaRNA Package\RNAfold.exe', 
                        help='RNAfold可执行文件路径')
    
    args = parser.parse_args()
    
    # 验证输入文件和RNAfold路径
    if not os.path.exists(args.input):
        print(f"输入文件不存在: {args.input}")
        return
    
    if not os.path.exists(args.rnafold):
        print(f"RNAfold可执行文件不存在: {args.rnafold}")
        return
    
    # 运行预测
    predict_secondary_structures(args.input, args.output, args.rnafold)

if __name__ == '__main__':
    main()