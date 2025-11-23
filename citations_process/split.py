import json
import re

# 输入和输出文件路径
input_json_file = 'F:\\citation\\extracted_data.json'  # 提取后的 JSON 文件
output_graph_file = 'F:\\citation\\graph.txt'          # 存储图结构的文件
output_attribute_file = 'F:\\citation\\user_attributes_int.txt'  # 存储论文属性的文件
output_attribute_inf_file = 'F:\\citation\\attribute.txt'  # 存储论文属性的文件
stem_path = "D:\\Desktop\\CS\\dataProcess\\stemmer.lowercase.txt"
stop_path = "D:\\Desktop\\CS\\dataProcess\\stopword.txt"
stem = {}
stop = set()
try:
    with open(stem_path, 'r') as f:
        for line in f:
            s = line.strip().split("/")
            for i in range(1, len(s)):
                stem[s[i]] = s[0]
    
except Exception as e:
    print(e)

try:
    with open(stop_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                stop.add(line)
    
except Exception as e:
    print(e)
# 读取提取的数据
with open(input_json_file, 'r', encoding='utf-8') as in_file:
    extracted_data = json.load(in_file)
def stemmer(word,stem):
    if word in stem:
        word = stem[word]
    return word
# 创建论文 ID 映射（原始 paper_id -> 新分配的连续 id）
paper_title_map = {}
for item in extracted_data:
    original_id = item['id']
    title = item['title']
    title = re.sub(r'[^a-zA-Z ]', ' ', title.lower())
    # title = ' '.join(title.split())  # 可选：压缩多余空格
    words = title.split()
    stemWord = set(stemmer(w,stem) for w in words if w)
    filtered = [w for w in stemWord if w not in stop]
    item['title'] = ' '.join(filtered)
    if re.search(r'[a-zA-Z]', item['title']):
        paper_title_map[original_id] = item['title']

    #     continue
    # if original_id not in paper_id_map:
    #     paper_id_map[original_id] = current_id
    #     current_id += 1
adj_list = {}
EN = 0
for item in extracted_data:
    if item['id'] not in paper_title_map:
        continue
    source_id = item['id']
    for ref_id in item.get('references', []):
        if ref_id in paper_title_map:
            target_id = ref_id

            # 添加双向边
            if source_id not in adj_list:
                adj_list[source_id] = set()
            adj_list[source_id].add(target_id)
            EN += 1

            if target_id not in adj_list:
                adj_list[target_id] = set()
            adj_list[target_id].add(source_id)
paper_id_map = {}

current_id = 0
# 收集所有非孤立节点
for node_id in adj_list:
    if len(adj_list[node_id]) > 0:
        paper_id_map[node_id] = current_id
        current_id += 1
attr_id_map = {}
current_id = 0
attr_length = 0
# 收集所有非孤立节点
for node_id in adj_list:
    if len(adj_list[node_id]) > 0:
        word = paper_title_map[node_id].split()
        attr_length += len(word)
        for w in word:
            if w not in attr_id_map:
                attr_id_map[w] = current_id
                current_id += 1

# 写入 attribute.txt 和 graph.txt
with open(output_attribute_inf_file, 'w', encoding='utf-8') as attr_inf_file, \
     open(output_attribute_file, 'w', encoding='utf-8') as attr_file, \
     open(output_graph_file, 'w', encoding='utf-8') as graph_file:
    graph_file.write(f"#Node Number: {len(paper_id_map)}\tEdge Number: {EN}\n")
    attr_file.write(f"#Node Number: {len(paper_id_map)}\tatt Number: {current_id}\tattLength: {attr_length}\n")
    # 遍历数据，写入属性和边
    written_edges = set()
    for src in adj_list:
        line = str(paper_id_map[src])  
        attid = list()
        for word in paper_title_map[src].split():
            line += '\t' + str(attr_id_map[word]) + "\t" + str(word)      
            attid.append(attr_id_map[word])    
        sorted_attid = sorted(attid)
        attr_inf_file.write(line + "\n")
        line2 = str(paper_id_map[src]) + '\t'
        for i in range(0, len(sorted_attid)-1):
            line2 += str(sorted_attid[i]) + ','
        line2 += str(sorted_attid[-1])
        attr_file.write(line2 + "\n")
        # attr_file.write(f"{paper_id_map[src]}\t{','.join(attid)}\n")
        for dst in adj_list[src]:
            edge = (src, dst)
            reverse_edge = (dst, src)
            if edge not in written_edges and reverse_edge not in written_edges:
                written_edges.add(edge)
                graph_file.write(f"{paper_id_map[src]}\t{paper_id_map[dst]}\n")

print(f"编号与属性已写入 {output_attribute_file}")
print(f"引用关系已写入 {output_graph_file}")