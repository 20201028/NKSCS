import random
import time
import os
import argparse

def generate_scaled_graph_simple(dataset: str, num_copies: int, node_num: int = 1000):
    """
    简化版的图扩展生成器
    """
    # 构建输入文件路径
    original_file = os.path.join("DataGraph", dataset, "user_attributes_int.txt")
    
    
    # 构建输出文件夹路径
    output_folder = os.path.join("DataGraph", f"{dataset}{num_copies}")
    output_file = os.path.join(output_folder, "user_attributes_int.txt")
    
    # 创建输出文件夹（如果不存在）
    os.makedirs(output_folder, exist_ok=True)
    
    print(f"生成 {num_copies} 倍扩展属性文件...")
    print(f"输入文件: {original_file}")
    print(f"输出文件: {output_file}")
    
    start_time = time.time()
    
    with open(output_file, 'w') as out_f:
        # 复制所有副本的内部属性
        for copy_id in range(num_copies):
            offset = copy_id * node_num
            print(f"处理副本 {copy_id + 1}/{num_copies}")
            
            with open(original_file, 'r') as in_f:
                for line in in_f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 2:
                            u = int(parts[0]) + offset
                            att = parts[1]
                            out_f.write(f"{u}\t{att}\n")
    
    generation_time = time.time() - start_time
    file_size = os.path.getsize(output_file) / (1024 * 1024 * 1024)
    
    print(f"生成完成!")
    print(f"文件大小: {file_size:.2f} GB")
    print(f"生成时间: {generation_time:.2f} 秒")

# 使用示例
if __name__ == "__main__":
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='生成扩展属性文件')
    parser.add_argument('--dataset', type=str, required=True, help='数据集名称')
    parser.add_argument('--scales', type=int, nargs='+', default=[20, 40, 60, 80, 100], help='扩展倍数列表')
    parser.add_argument('--node_num', type=int, default=1278674, help='每个副本的节点数')
    
    args = parser.parse_args()
    
    # 生成多个规模的属性文件
    for scale in args.scales:
        generate_scaled_graph_simple(
            dataset=args.dataset,
            num_copies=scale,
            node_num=args.node_num
        )