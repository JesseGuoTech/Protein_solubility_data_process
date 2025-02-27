from Bio import SeqIO
import random


def split_fasta(input_file, output_files, proportions):
    """
    将FASTA文件按比例拆分为多个子文件。

    参数：
    - input_file: 输入FASTA文件路径
    - output_files: 输出文件列表（如 ['train.fasta', 'test.fasta', 'nesg.fasta', 'chang.fasta']）
    - proportions: 对应的比例列表（如 [0.7, 0.15, 0.075, 0.075]）
    """
    # 读取FASTA文件中的记录
    records = list(SeqIO.parse(input_file, "fasta"))

    # 打乱顺序
    random.shuffle(records)

    # 计算各部分数量
    total = len(records)
    sizes = [int(total * prop) for prop in proportions]

    # 修正最后一部分的大小以确保总数一致
    sizes[-1] = total - sum(sizes[:-1])

    # 按大小切分记录
    splits = []
    start = 0
    for size in sizes:
        splits.append(records[start:start + size])
        start += size

    # 将切分的记录写入对应文件
    for output_file, split in zip(output_files, splits):
        SeqIO.write(split, output_file, "fasta")
    print("FASTA文件已成功拆分！,拆分为", sizes)


# 输入和输出文件路径
input_fasta = "train_origin.fasta"  # 你的长FASTA文件
output_files = ["train.fasta", "test.fasta", "val1.fasta", "val2.fasta"]
proportions = [0.7, 0.15, 0.075, 0.075]

# 拆分FASTA文件
split_fasta(input_fasta, output_files, proportions)
