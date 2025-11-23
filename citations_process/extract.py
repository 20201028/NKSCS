import json

# 输入和输出文件路径
input_json_file = 'F:\citation\dblp_v14.json'  # 替换为你的 JSON 文件路径
output_json_file = 'F:\citation\extracted_data.json'  # 提取的 JSON 数据将被写入这个文件

# 打开 JSON 文件并加载数据
with open(input_json_file, 'r', encoding='utf-8') as json_file:
    data = json.load(json_file)

# 用于存储提取的数据
extracted_data = []

# 遍历 JSON 数组中的每个对象
for item in data:
    paper_id = item.get('id', 'N/A')
    title = item.get('title', 'N/A')
    references = item.get('references', [])

    # 构建新的 JSON 对象
    extracted_item = {
        "id": paper_id,
        "title": title,
        "references": references
    }

    extracted_data.append(extracted_item)

# 将提取的数据写入新的 JSON 文件
with open(output_json_file, 'w', encoding='utf-8') as out_file:
    json.dump(extracted_data, out_file, indent=4)

print(f"提取完成，结果已写入 {output_json_file}")