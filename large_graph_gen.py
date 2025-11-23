import random
import time
import os
import argparse

def generate_scaled_graph_simple(dataset: str, num_copies: int, connection_prob: float = 0.001, node_num: int = 1000):
    """
    简化版的图扩展生成器
    """
    # 构建输入文件路径
    original_file = os.path.join("DataGraph", dataset, "graph.txt")
    
    # 构建输出文件夹路径
    output_folder = os.path.join("DataGraph", f"{dataset}{num_copies}")
    output_file = os.path.join(output_folder, "graph.txt")
    
    # 创建输出文件夹（如果不存在）
    os.makedirs(output_folder, exist_ok=True)
    
    print(f"生成 {num_copies} 倍扩展图...")
    print(f"输入文件: {original_file}")
    print(f"输出文件: {output_file}")
    
    start_time = time.time()
    connections_made = 0
    
    with open(output_file, 'w') as out_f:
        # 复制所有副本的内部边
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
                            v = int(parts[1]) + offset
                            out_f.write(f"{u}\t{v}\n")
        
        # 添加连接边
        print("添加连接边...")
        for i in range(num_copies):
            for j in range(i + 1, num_copies):
                offset_i = i * node_num
                offset_j = j * node_num
                
                # 修正连接边生成逻辑
                for _ in range(int(node_num * connection_prob)):
                    node_i = random.randint(0, node_num-1) + offset_i
                    node_j = random.randint(0, node_num-1) + offset_j
                    out_f.write(f"{node_i}\t{node_j}\n")
                    connections_made += 1
    
    generation_time = time.time() - start_time
    file_size = os.path.getsize(output_file) / (1024 * 1024 * 1024)
    
    print(f"生成完成!")
    print(f"连接边数量: {connections_made}")
    print(f"文件大小: {file_size:.2f} GB")
    print(f"生成时间: {generation_time:.2f} 秒")

# 使用示例
if __name__ == "__main__":
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='生成扩展图')
    parser.add_argument('--dataset', type=str, required=True, help='数据集名称')
    parser.add_argument('--scales', type=int, nargs='+', default=[20, 40, 60, 80, 100], help='扩展倍数列表')
    parser.add_argument('--connection_prob', type=float, default=0.001, help='连接概率')
    parser.add_argument('--node_num', type=int, default=1278674, help='每个副本的节点数')
    
    args = parser.parse_args()
    
    # 生成多个规模的图
    for scale in args.scales:
        generate_scaled_graph_simple(
            dataset=args.dataset,
            num_copies=scale,
            connection_prob=args.connection_prob,
            node_num=args.node_num
        )