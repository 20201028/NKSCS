import xml.etree.ElementTree as ET
# ext_path = "D:\\Desktop\\CS\\dataProcess\\dblp-2024-06-02-ext.xml"
# txt_path = "D:\\Desktop\\CS\\dataProcess\\dblp-2024-06-02-txt.txt"
ext_path = "D:\\Desktop\\CS\\dataProcess\\dblp-ext.xml"
txt_path = "D:\\Desktop\\CS\\dataProcess\\dblp-txt.txt"
    
try:
    # 解析XML文件
    tree = ET.parse(ext_path)
    root = tree.getroot()
        
    print(f"nodeList.getLength(): {len(root)}")
        
    # 打开输出文件进行写入
    with open(txt_path, 'w', encoding='ISO-8859-1') as stdout:
        # 遍历XML的根元素的子节点
        for node in root:
            node_name = node.tag

            title = ""
            author_list = []

            title_element = node.find('title')

            # Extract the text content including nested tags
            if title_element is not None:
                title = ''.join(title_element.itertext())
            
                    
            # 遍历文章节点的子节点，获取标题和作者

            for child in node:
                # if 'title' == child.tag:
                #     title = child.text.strip()
                if 'author' == child.tag:
                    author_list.append(child.text.strip())
                    
            # 输出标题和作者信息（跳过没有作者的文章）
            if author_list:
                stdout.write(title + '\n')
                stdout.write(author_list[0])
                for author in author_list[1:]:
                    stdout.write('\t' + author)
                stdout.write('\n')

except Exception as e:
    print(f"An error occurred: {e}")