import sys
import os
import glob
import subprocess
import gzip
import shutil

def run_cmd(cmd):
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def prepare_fastq_files(input_dir):
    """
    Quét thư mục, tìm file .fq, .fastq, .fq.gz, .fastq.gz.
    Nếu là file nén .gz, tiến hành giải nén. Trả về danh sách file .fastq.
    """
    ready_fastq_files = []
    
    # Tìm tất cả các file trong thư mục input
    all_files = glob.glob(os.path.join(input_dir, "*"))
    
    for file_path in all_files:
        file_name = os.path.basename(file_path).lower()
        
        # 1. Nếu là file nén (.fq.gz hoặc .fastq.gz)
        if file_name.endswith('.fq.gz') or file_name.endswith('.fastq.gz'):
            # Tạo tên file mới bỏ đuôi .gz và chuẩn hóa thành .fastq
            new_name = file_name.replace('.fq.gz', '.fastq').replace('.fastq.gz', '.fastq')
            uncompressed_path = os.path.join(input_dir, new_name)
            
            # Chỉ giải nén nếu file .fastq chưa tồn tại (tiết kiệm thời gian chạy lại)
            if not os.path.exists(uncompressed_path):
                print(f"📦 Đang giải nén file: {file_name} -> {new_name}")
                with gzip.open(file_path, 'rb') as f_in:
                    with open(uncompressed_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
            else:
                print(f"✅ Đã có sẵn file giải nén cho: {file_name}")
                
            ready_fastq_files.append(uncompressed_path)
            
        # 2. Nếu đã là file chưa nén (.fq hoặc .fastq)
        elif file_name.endswith('.fastq') or file_name.endswith('.fq'):
            ready_fastq_files.append(file_path)

    return ready_fastq_files

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python pipeline.py <input_dir> <output_dir>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    # Tạo các thư mục con (B1, B2, B3)
    qc_dir = os.path.join(output_dir, "B1_NanoPlot")
    filt_dir = os.path.join(output_dir, "B2_NanoFilt")
    emu_dir = os.path.join(output_dir, "B3_EMU")
    os.makedirs(qc_dir, exist_ok=True)
    os.makedirs(filt_dir, exist_ok=True)
    os.makedirs(emu_dir, exist_ok=True)

    print("🔍 Đang kiểm tra và chuẩn bị dữ liệu đầu vào...")
    # Tự động nhận diện và giải nén
    fastq_files = prepare_fastq_files(input_dir)
    
    if not fastq_files:
        print("❌ Lỗi: Không tìm thấy file dữ liệu nào (fq, fastq, fq.gz) trong thư mục input!")
        sys.exit(1)

    print(f"🚀 Bắt đầu xử lý {len(fastq_files)} mẫu dữ liệu...")

    # =========================================================================
    # Bước 1 & Bước 2 & Bước 3.1: NẰM TRONG VÒNG LẶP (Chạy riêng cho TỪNG mẫu)
    # =========================================================================
    for fq in fastq_files:
        # Xóa các đuôi mở rộng để lấy tên Mẫu chuẩn (SampleID)
        sample_name = os.path.basename(fq)
        for ext in ['.fastq', '.fq', '.fastq.gz', '.fq.gz']:
            sample_name = sample_name.replace(ext, "")
            
        filtered_fq = os.path.join(filt_dir, f"{sample_name}_filtered.fastq")
        
        # B1: NanoPlot
        run_cmd(f"NanoPlot -t 4 --fastq {fq} -o {qc_dir}/{sample_name}")
        
        # B2: NanoFilt (Sử dụng luồng Python an toàn thay cho dấu < > của Shell)
        print(f"Running: NanoFilt cho mẫu {sample_name}...")
        with open(fq, "r") as f_in, open(filtered_fq, "w") as f_out:
            subprocess.run(
                ["NanoFilt", "-q", "10", "-l", "1400", "--maxlength", "1700"],
                stdin=f_in, stdout=f_out, check=True
            )
        
        # B3.1: EMU Abundance
        run_cmd(f"emu abundance {filtered_fq} --db emu_db_mar2026 --output-dir {emu_dir}/{sample_name} --threads 4")

    # =========================================================================
    # CÁC BƯỚC DƯỚI ĐÂY NẰM NGOÀI VÒNG LẶP (Chỉ chạy 1 lần sau khi gộp đủ 36 mẫu)
    # =========================================================================
    
    print("🔄 Đang gộp kết quả EMU của tất cả các mẫu...")
    # B3.2: Gộp kết quả EMU thành ma trận tổng
    run_cmd(f"emu combine --dir {emu_dir} --output {emu_dir}/combined_abundance_matrix.tsv")

    print("📊 Đang lọc ma trận (>= 2 samples, >= 0.4%)...")
    # B3.3: Gọi Python lọc điều kiện (Sử dụng sys.executable để đảm bảo gọi đúng môi trường)
    raw_matrix = os.path.join(emu_dir, "combined_abundance_matrix.tsv")
    filtered_matrix = os.path.join(output_dir, "filtered_matrix.csv")
    run_cmd(f"{sys.executable} filter_matrix.py {raw_matrix} {filtered_matrix}")

    print("📈 Đang chạy phân tích thống kê sinh thái bằng R...")
    # B4: Gọi script R thống kê sinh thái
    metadata_file = os.path.join(input_dir, "metadata.csv")
    run_cmd(f"Rscript analysis.R {filtered_matrix} {metadata_file} {output_dir}")

    print("🎉 Pipeline completed successfully!")