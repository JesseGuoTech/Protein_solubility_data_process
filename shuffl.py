from Bio import SeqIO
import random


def split_fasta(input_file, train_file, validation_file, validation_percentage=0.03):
    """
    读取FASTA文件，打乱顺序并拆分为训练集和验证集。

    参数：
    - input_file: 输入的FASTA文件路径
    - train_file: 训练集FASTA文件路径
    - validation_file: 验证集FASTA文件路径
    - validation_percentage: 验证集比例（默认30%）
    """
    # 读取FASTA文件中的所有记录
    records = list(SeqIO.parse(input_file, "fasta"))

    # 打乱顺序
    random.shuffle(records)

    # 计算验证集数量
    validation_count = int(len(records) * validation_percentage)

    # 拆分数据集
    validation_records = records[:validation_count]  # 前30%
    train_records = records[validation_count:]  # 剩余70%

    # 保存训练集和验证集
    SeqIO.write(train_records, train_file, "fasta")
    SeqIO.write(validation_records, validation_file, "fasta")

    print(f"处理完成！")
    print(f"训练集保存为: {train_file}，包含 {len(train_records)} 条记录")
    print(f"验证集保存为: {validation_file}，包含 {len(validation_records)} 条记录")


# 输入和输出文件路径
input_file = "all.fasta"
train_file = "train_origin.fasta"
validation_file = "validation.fasta"

# 调用函数进行处理
split_fasta(input_file, train_file, validation_file)
