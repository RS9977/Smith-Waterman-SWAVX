def are_files_same(file1_path, file2_path):
    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
        for line1, line2 in zip(file1, file2):
            if line1.strip() != line2.strip():
                return False

        # Check if the files have different number of lines
        remaining_lines_file1 = list(file1)
        remaining_lines_file2 = list(file2)
        if remaining_lines_file1 or remaining_lines_file2:
            return False

    return True

if __name__ == "__main__":
    file1_path = "SWalgo.txt"
    file2_path = "256_ST_SB.txt"

    if are_files_same(file1_path, file2_path):
        print("The files are the same.")
    else:
        print("The files are different.")
