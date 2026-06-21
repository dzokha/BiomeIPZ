# BiomeIPZ - Automated Nanopore 16S rRNA Microbiome Pipeline

**BiomeIPZ** là một hệ thống pipeline tin sinh học tự động hóa, được thiết kế chuyên dụng để xử lý dữ liệu giải trình tự gen thế hệ mới (NGS) từ nền tảng **Oxford Nanopore Technologies (ONT)** nhắm vào vùng gen **16S rRNA**. Hệ thống tích hợp một giao diện Web API (FastAPI) giúp quản lý và kích hoạt các tiến trình chạy ngầm, từ bước kiểm định chất lượng thô cho đến phân tích thống kê sinh thái.

---

## 💻 1. Môi trường Hệ thống (System Environment)

Pipeline được tối ưu hóa tốt nhất trên hệ điều hành **macOS** sử dụng môi trường ảo cô lập **Miniconda** để tránh xung đột thư viện.

* **Python Version:** `3.11` (Conda environment: `biome_env`)
* **Các thư viện lõi:** `fastapi`, `uvicorn`, `pydantic`, `pandas`, `matplotlib (< 3.9.0)`, `ete3`
* **Công cụ Tin sinh học:** `NanoPlot`, `NanoFilt`, `EMU`, `R-base`

---

## 🏗️ 2. Kiến trúc Pipeline (Pipeline Architecture)

Pipeline vận hành tuần tự qua 4 bước cốt lõi:
1.  **B1_NanoPlot:** Kiểm tra và trực quan hóa chất lượng của các luồng đọc (reads) thô từ file `.fastq` hoặc `.fq.gz`.
2.  **B2_NanoFilt:** Cắt lọc dữ liệu thô bằng luồng native Python. Loại bỏ reads có chất lượng thấp (Q < 10) và giới hạn độ dài reads tối ưu cho vùng 16S (từ 1400 bp đến 1700 bp).
3.  **B3_EMU (Taxonomic Profiling):** Định danh loài vi khuẩn bằng thuật toán EMU Abundance, sau đó gộp kết quả các mẫu thành một ma trận tổng (`combined_abundance_matrix.tsv`).
4.  **Statistical Analysis:** Chạy script Python (`filter_matrix.py`) để lọc nhiễu (xuất hiện ở >= 2 samples và đạt tỉ lệ >= 0.4%) và gọi `Rscript analysis.R` để vẽ biểu đồ thống kê sinh thái.

---

## 🗄️ 3. Cẩm nang Xây dựng Cơ sở dữ liệu EMU từ con số 0
(Custom EMU Database Generation Guide)

EMU yêu cầu một bộ cơ sở dữ liệu (Database) tham chiếu có cấu trúc khắt khe bao gồm file trình tự FASTA và file phả hệ Taxonomy tương ứng. Quy trình dưới đây giúp xây dựng một Database tối tân nhất (Cập nhật phiên bản mới nhất năm 2026).

### Nguyên liệu thô (Tải thủ công từ NCBI)
1. **File trình tự 16S RefSeq Nucleotide (sequence.fasta — 41.9 MB)**
- Mục đích: Chứa bản đồ ADN thô (chuỗi A, T, G, C) của toàn bộ 26.870 chủng vi khuẩn trên thế giới để làm dữ liệu đối chiếu.
- Các bước tải: 
+ Bước 1: Truy cập vào trang quản lý của NCBI qua đường link: https://www.ncbi.nlm.nih.gov/refseq/targetedloci/16S_process/
+ Bước 2: Click chuột vào dòng chữ "16S RefSeq Nucleotide sequence records".
+ Bước 3: Tại trang danh sách hiện ra, bạn nhìn sang góc trên cùng bên phải, click vào nút Send to $\rightarrow$ chọn ô File $\rightarrow$ tại mục Format chọn FASTA $\rightarrow$ bấm nút Create File để tải về một file duy nhất.
2. **File Danh bạ Accession Number (RefSeq-release235.catalog.gz — 3.3 GB)**
- Mục đích: Đóng vai trò như một cuốn "từ điển" ngoại tuyến, giúp phần mềm tự động dịch từ mã quản lý khoa học của NCBI (ví dụ: NR_113675.1) sang mã số định danh loài (TaxID) để tra cứu phả hệ.
- Các bước tải:
+ Bước 1: Truy cập vào kho máy chủ FTP của NCBI qua đường link: https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/
+ Bước 2: Cuộn chuột tìm đúng file có tên chính xác là RefSeq-release235.catalog.gz.
- Bước 3: Click chuột thẳng vào tên file đó để tải về máy tính (Lưu ý: Giữ nguyên file nén đuôi .gz, tuyệt đối không giải nén).

### Quy trình xử lý Offline (Offline Mapping Workflow)

Do NCBI chặn các truy vấn online số lượng lớn, quy trình được thực hiện hoàn toàn ngoại tuyến thông qua script `build_emu_raw.py` bằng cơ chế đọc luồng nén (Stream decompression) để tiết kiệm tài nguyên RAM.

#### Bước 3.1: Tạo file dữ liệu sạch bằng Script Python (`build_emu_raw.py`)

1. **Chuẩn bị thư mục:** Bạn tạo một thư mục mới tên là `db` nằm ngay bên trong thư mục chứa file script `build_emu_raw.py`.
2. **Đặt file vào vị trí:** Di chuyển (hoặc copy) cả 2 file nguyên liệu là `sequence.fasta` và `RefSeq-release235.catalog.gz` bỏ vào bên trong thư mục `db` vừa tạo ở trên.
3. **Khởi chạy lệnh:** Bạn mở Terminal, di chuyển vào thư mục chứa file script và gõ câu lệnh sau để bắt đầu tiến trình xử lý:

```bash
python build_emu_raw.py
```
Cơ chế: Script tự động quét file FASTA để bóc tách mã Accession, sau đó dò tìm trong file Catalog dung lượng 3.3 GB để lấy TaxID tương ứng. Cuối cùng, thư viện ete3 sẽ kết nối với cây taxonomy cục bộ để dịch ngược ra đủ 7 cấp bậc phân loại sinh giới.

Kết quả: Sinh ra 2 file chuẩn chỉnh: species_taxid.fasta (Trình tự gen sạch) và taxonomy.tsv (Sơ đồ phả hệ sinh giới).

#### Bước 3.2: Lập chỉ mục bằng EMU (Database Indexing)
Tạo một thư mục mới tên là my_emu_db_2026, chuyển 2 file vừa sinh ra vào thư mục này và chạy lệnh đóng gói của EMU:

```bash
emu build-database --dir my_emu_db_2026
```
Hệ thống sẽ kích hoạt minimap2 ngầm để chia nhỏ trình tự và tạo ma trận tìm kiếm siêu tốc (file .mmi).

## 🚀 4. Hướng dẫn Khởi chạy Dự án (How to Run)

### Cách 1: Chạy qua Giao diện Web API (Swagger UI)

1. Kích hoạt môi trường và bật server FastAPI:

```bash
conda activate biome_env
cd /Users/dzokha/Desktop/Project/BiomeIPZ
uvicorn main:app --host 0.0.0.0 --port 8000
```

2. Truy cập vào trình duyệt: http://127.0.0.1:8000/docs

3. Tìm đến API /run-pipeline, bấm Try it out và điền dữ liệu JSON với đường dẫn tuyệt đối chuẩn macOS (Sử dụng dấu gạch chéo xuôi /):

```json
{
  "job_id": "TEST_36_SAMPLES",
  "input_folder": "/Users/dzokha/Desktop/Project/BiomeIPZ/data/test_run_01",
  "output_folder": "/Users/dzokha/Desktop/Project/BiomeIPZ/data/results_01"
}
```
4. Bấm Execute. Tiến trình sẽ tự động chạy ngầm (Background Task). Bạn có thể theo dõi tiến độ thời gian thực bằng cách đọc file log tại thư mục kết quả.

### Cách 2: Chạy trực tiếp qua Terminal (CLI Mode)
Nếu muốn ép lộ lỗi ngay trên màn hình để kiểm tra mã nguồn, chạy lệnh trực tiếp:

```bash
python3 pipeline.py /Users/dzokha/Desktop/Project/BiomeIPZ/data/test_run_01 /Users/dzokha/Desktop/Project/BiomeIPZ/data/results_01
```

## 📂 5. Cấu trúc Thư mục Dự án (Directory Structure)

```text
BiomeIPZ/
├── main.py                     # Web API Server (FastAPI)
├── pipeline.py                 # Script điều tốc và kích hoạt Pipeline chính
├── filter_matrix.py            # Script Python lọc ma trận kết quả
├── analysis.R                  # Script R chạy thống kê sinh thái và phân tích biểu đồ
├── build_emu_raw.py            # Script tự động đúc Custom Database ngoại tuyến
├── requirements.txt            # Danh sách các thư viện Python cần thiết
├── my_emu_db_2026/             # Thư mục Cơ sở dữ liệu EMU sau khi build thành công
│   ├── species_taxid.fasta
│   └── taxonomy.tsv
└── data/
    ├── test_run_01/            # Thư mục dữ liệu đầu vào (Chứa file .fq.gz & metadata.csv)
    └── results_01/             # Thư mục đầu ra chứa kết quả tự động phân tách B1, B2, B3
```