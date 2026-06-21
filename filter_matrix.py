import pandas as pd
import sys

if len(sys.argv) != 3:
    print("Usage: python filter_matrix.py <input_emu_tsv> <output_filtered_csv>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Đọc ma trận từ EMU (thường EMU có cột 'tax_id', 'abundance', v.v.)
# Giả định ma trận định dạng: Hàng = Loài (Taxa), Cột = Tên các mẫu
df = pd.read_csv(input_file, sep='\t', index_col=0)

# Lọc bỏ các cột metadata thừa của công cụ nếu có (chỉ giữ lại các cột chứa số liệu)
numeric_df = df.select_dtypes(include='number')

# Tính toán Độ phong phú tương đối (Relative Abundance) cho toàn bộ ma trận
relative_abundance_df = numeric_df.div(numeric_df.sum(axis=0), axis=1)

# Điều kiện 1: Xuất hiện ở >= 2 mẫu (Giá trị > 0)
mask_samples = (numeric_df > 0).sum(axis=1) >= 2

# Điều kiện 2: Độ phong phú cực đại >= 0.4% (tức là 0.004)
mask_abundance = relative_abundance_df.max(axis=1) >= 0.004

# Áp dụng cả 2 điều kiện (AND)
final_mask = mask_samples & mask_abundance
filtered_df = df[final_mask]

# Xuất kết quả
filtered_df.to_csv(output_file)
print(f"Lọc thành công: Từ {len(df)} loài xuống còn {len(filtered_df)} loài.")