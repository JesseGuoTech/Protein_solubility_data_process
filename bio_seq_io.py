from pathlib import Path
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
import csv

csv_files_path = "csv"
fasta_files_path = "fasta"
out_files_dir = "output"
sol_border: float = 0  # 是否可溶的临界值
files_directory: Path = Path(__file__).parent

def aa_replace(csv_files_p, fasta_files_p, out_files_p, border) -> tuple:
    """
    read csv mutant files and original protein fasta files , convert them into fasta files with solubility(0,1)
    Args:
        csv_files_p:csv mutant input file path
        fasta_files_p: fasta input files path
        out_files_p: output file path
        border: The ddg border to judge if protein is soluble

    Returns:
        A tuple with two elements:File with errors(return None if no error) & processing status
    """
    csv_files, fasta_files, out_files = files_directory.joinpath(csv_files_p), files_directory.joinpath(
        fasta_files_p), files_directory.joinpath(out_files_p)
    out_files.mkdir(parents=True, exist_ok=True)

    for csv_file in csv_files.glob("*.csv"):
        out_file = out_files.joinpath(csv_file.stem + ".fasta")
        with open(csv_file, "r") as csv_obj:
            # 打开文件 读取csv 读取序列
            csv_reader = csv.reader(csv_obj)
            seq = SeqIO.read(fasta_files.joinpath(csv_file.with_suffix(".fasta").name), "fasta").seq
            mod_seq = []
            for serial_number, resident, origin_aa, turned_aa, _, _, ddg in csv_reader:
                sol = 1 if float(ddg) < border else 0
                mul_seq = MutableSeq(seq)
                if mul_seq[int(resident) - 1] == origin_aa:
                    mul_seq[int(resident) - 1] = turned_aa
                else:
                    return out_file.name, False
                mod_seq.append(SeqRecord(seq=mul_seq, id=serial_number, description="ddg=" + ddg + "|" + str(sol)))
            SeqIO.write(mod_seq, out_file, "fasta")
    return None, True


if __name__ == "__main__":
    fail_file_name, run_statue = aa_replace(csv_files_path, fasta_files_path, out_files_dir, sol_border)
    if run_statue:
        print("Done")
    else:
        print("Error,file name = ", fail_file_name)
