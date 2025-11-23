import json

# 输入和输出文件路径
input_json_file = 'F:\\citation\\extracted_data.json'  # 已提取的 JSON 文件
output_json_file = 'F:\\citation\\filtered_data.json'  # 过滤后的输出文件

# 读取提取的 JSON 数据
with open(input_json_file, 'r', encoding='utf-8') as in_file:
    data = json.load(in_file)

# 过滤掉 references 为空的条目
filtered_data = [item for item in data if item.get("references")]

# 将过滤后的数据写入新的 JSON 文件
with open(output_json_file, 'w', encoding='utf-8') as out_file:
    json.dump(filtered_data, out_file, indent=4)

print(f"过滤完成，结果已写入 {output_json_file}")