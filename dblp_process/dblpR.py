import xml.etree.ElementTree as ET

def main():
    # 定义输入和输出文件路径
    # in_path = "D:\\Desktop\\CS\\dataProcess\\dblp-2024-06-02.xml"
    # out_path = "D:\\Desktop\\CS\\dataProcess\\dblp-2024-06-02-R.xml"
    in_path = "D:\\Desktop\\CS\\dataProcess\\dblp.xml"
    out_path = "D:\\Desktop\\CS\\dataProcess\\dblp-R.xml"

    try:
        # Open input file for reading
        with open(in_path, 'r', encoding='ISO-8859-1') as stdin:
            # Open output file for writing
            with open(out_path, 'w', encoding='ISO-8859-1') as stdout:
                for line in stdin:
                    if('&' in line):
                        line = line.replace('&', '&amp;')
                    # if('&uuml;' in line):
                    #     line = line.replace('&uuml;', 'ü')
                    # if('&uacute;' in line):
                    #     line = line.replace('&uacute;', 'ú')
                    # if('&auml;' in line):
                    #     line = line.replace('&auml;', 'ä')
                    # if('&aacute;' in line):
                    #     line = line.replace('&aacute;', 'á')
                    # if('&atilde;' in line):
                    #     line = line.replace('&atilde;', 'ã')
                    # if('&aelig;' in line):
                    #     line = line.replace('&aelig;', 'æ')
                    # if('&Aacute;' in line):
                    #     line = line.replace('&Aacute;', 'Á')
                    # if('&ouml;' in line):
                    #     line = line.replace('&ouml;', 'ö')
                    # if('&oacute;' in line):
                    #     line = line.replace('&oacute;', 'ó')
                    # if('&ccedil;' in line):
                    #     line = line.replace('&ccedil;', 'ç')
                    # if('&Ccedil;' in line):
                    #     line = line.replace('&Ccedil;', 'Ç')
                    # if('&iacute;' in line):
                    #     line = line.replace('&iacute;', 'í')
                    # if('&icirc;' in line):
                    #     line = line.replace('&icirc;', 'î')
                    # if('&iuml;' in line):
                    #     line = line.replace('&iuml;', 'ï')
                    # if('&egrave;' in line):
                    #     line = line.replace('&egrave;', 'è')
                    # if('&eacute;' in line):
                    #     line = line.replace('&eacute;', 'é')
                    # if('&euml;' in line):
                    #     line = line.replace('&euml;', 'ë')
                    # if('&reg;' in line):
                    #     line = line.replace('&reg;', '®')
                    # if('&ntilde;' in line):
                    #     line = line.replace('&ntilde;', 'ñ')
                    # if('&oslash;' in line):
                    #     line = line.replace('&oslash;', 'ø')
                    # if('&szlig;' in line):
                    #     line = line.replace('&szlig;', 'ß')
                    stdout.write(line)
    except IOError as e:
        print(f"文件操作错误: {e}")

if __name__ == "__main__":
    main()