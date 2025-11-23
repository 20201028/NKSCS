import xml.etree.ElementTree as ET

def main():
    # 定义输入和输出文件路径
    # in_path = "D:\\Desktop\\CS\\dataProcess\\dblp-2024-06-02-R.xml"
    # out_path = "D:\\Desktop\\CS\\dataProcess\\dblp-2024-06-02-ext.xml"
    in_path = "D:\\Desktop\\CS\\dataProcess\\dblp-R.xml"
    out_path = "D:\\Desktop\\CS\\dataProcess\\dblp-ext.xml"

    try:
        # 解析XML文件
        tree = ET.parse(in_path)
        root = tree.getroot()

        with open(out_path, 'w', encoding='ISO-8859-1') as stdout:
            # 写入XML头
            stdout.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
            stdout.write('<!DOCTYPE dblp SYSTEM "dblp.dtd">\n')
            stdout.write('<dblp>\n')

            # 遍历根元素的子节点并筛选指定标签
            for element in root:
                tag = element.tag
                if tag in ["article", "inproceedings"]:
                # if tag in ["article"]:
                    stdout.write(ET.tostring(element, encoding='unicode'))

            # 写入XML尾部
            stdout.write('</dblp>\n')

    except ET.ParseError as e:
        print(f"XML解析错误: {e}")
    except IOError as e:
        print(f"文件操作错误: {e}")

if __name__ == "__main__":
    main()