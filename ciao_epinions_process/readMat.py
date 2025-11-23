import numpy as np
from scipy.io import loadmat

# 读取.mat文件
# data = loadmat('epinions\\rating.mat')
data = loadmat('ciao\\rating.mat')

# 提取rating数据
rating = data['rating']

# 初始化集合来存储唯一的用户和项
users = set()
items = set()
user_items = {}
# 输出文件路径
output_file = 'user_item.txt'
count = 0
# 将用户及其购买项写入txt文件，并统计用户数和项数
with open(output_file, 'w') as f:
    for row in rating:
        user_id, item_id = row[0], row[1]
        if user_id not in user_items:
            user_items[user_id] = set()
        user_items[user_id].add(item_id)
        users.add(user_id)
        items.add(item_id)
    
    # 统计用户个数和项个数
    userNum = len(users)
    itemNum = len(items)
    
    # 将统计信息写入txt文件开头
    f.write(f"#userNum: {userNum}, itemNum: {itemNum}\n")
    
   # 将用户及其购买项写入txt文件
    lines = []
    for user_id, item_id in user_items.items():
        count+=1
        if count == 40:
            sta= 0
        list_items = list(item_id)
        list_items.sort()
        # 确保 item_id 中的元素是字符串类型
        item_str = ",".join(map(str, list_items))
        lines.append(f"{user_id}\t{item_str}\n")
    
    # 批量写入文件
    f.writelines(lines)

# 读取.mat文件
# data = loadmat('epinions\\trustnetwork.mat')
data = loadmat('ciao\\trustnetwork.mat')

# 提取trustnetwork数据
trustnetwork = data['trustnetwork']

# 初始化集合来存储唯一的节点和边
users = set()
edges = set()

# 遍历trustnetwork数据
for row in trustnetwork:
    user1, user2 = row[0], row[1]
    users.add(user1)
    users.add(user2)
    # 将边存储为一个元组，确保无向边的顺序一致
    edge = tuple(sorted((user1, user2)))
    edges.add(edge)

# 统计节点个数和边个数
userNum = len(users)
edgeNum = len(edges)

# 输出文件路径
output_file = 'graph.txt'

# 将统计信息和边写入txt文件
with open(output_file, 'w') as f:
    f.write(f"#userNum: {userNum}, edgeNum: {edgeNum}\n")
    for edge in edges:
        f.write('\t'.join(map(str, edge)) + '\n')