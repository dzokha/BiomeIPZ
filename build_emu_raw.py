import os
import sys
import json
import urllib.request
import time
from ete3 import NCBITaxa

BATCH_SIZE = 300  # Gom cụm 300 mã gửi lên NCBI 1 lần để chạy siêu tốc

def fetch_taxids_from_ncbi(accessions):
    """ Gọi API NCBI E-utilities để đổi danh sách Accession -> TaxID """
    # Loại bỏ phiên bản nếu có (ví dụ: NR_113675.1 -> NR_113675) để tra cứu chính xác
    clean_accs = [acc.split('.')[0] for acc in accessions]
    ids_str = ",".join(clean_accs)
    
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={ids_str}&retmode=json"
    mapping = {}
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=20) as response:
            data = json.loads(response.read().decode())
            results = data.get("result", {})
            for uid in results.get("uids", []):
                uid_info = results.get(uid, {})
                acc_ver = uid_info.get("accessionversion")
                caption = uid_info.get("caption")
                taxid = uid_info.get("taxid")
                
                if taxid:
                    if acc_ver: mapping[acc_ver] = int(taxid)
                    if caption: mapping[caption] = int(taxid)
    except Exception as e:
        print(f"\n⚠️ Cảnh báo lỗi API NCBI: {e}. Sẽ thử lại...")
        time.sleep(1)
    return mapping

def build_emu_database_smart(input_fasta, out_fasta, out_tsv):
    print("⏳ Bước 1: Quét file FASTA để thu thập tất cả mã Accession...")
    accessions = []
    with open(input_fasta, 'r') as f:
        for line in f:
            if line.startswith(">"):
                # Lấy chữ đầu tiên sau dấu >, ví dụ: "NR_113675.1"
                acc = line[1:].split()[0]
                accessions.append(acc)
                
    total_records = len(accessions)
    print(f"✅ Tìm thấy {total_records} trình tự gen trong file thô.")

    print("\n🌐 Bước 2: Tự động kết nối NCBI lấy mã TaxID (Chạy ngầm theo block)...")
    acc_to_taxid = {}
    for i in range(0, total_records, BATCH_SIZE):
        batch = accessions[i:i+BATCH_SIZE]
        sys.stdout.write(f"\r⏳ Đang tải dữ liệu hộ chiếu sinh học: {i}/{total_records} con...")
        sys.stdout.flush()
        
        # Thử lại tối đa 3 lần nếu mạng lỗi
        for _ in range(3):
            mapping = fetch_taxids_from_ncbi(batch)
            if mapping:
                acc_to_taxid.update(mapping)
                break
        time.sleep(0.3) # Nghỉ một chút chống bị NCBI ban IP
    print(f"\n✅ Đã đồng bộ thành công TaxID cho các chủng vi khuẩn!")

    print("\n🌳 Bước 3: Khởi tạo dữ liệu cây phả hệ sinh giới...")
    ncbi = NCBITaxa()
    
    print("🧬 Bước 4: Đúc file species_taxid.fasta và taxonomy.tsv đúng chuẩn EMU...")
    with open(input_fasta, 'r') as f_in, \
         open(out_fasta, 'w') as f_fasta, \
         open(out_tsv, 'w') as f_tsv:
        
        # Ghi Header chuẩn của EMU
        f_tsv.write("tax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tclade\n")
        
        processed_taxids = set()
        current_acc = None
        valid_count = 0

        for line in f_in:
            line = line.strip()
            if line.startswith(">"):
                current_acc = line[1:].split()[0]
                
                # Tìm TaxID tương ứng
                tax_id_int = acc_to_taxid.get(current_acc)
                if not tax_id_int:
                    # Thử tìm bằng mã không chứa đuôi .1
                    tax_id_int = acc_to_taxid.get(current_acc.split('.')[0])

                if not tax_id_int:
                    current_acc = None  # Đánh dấu bị lỗi để bỏ qua chuỗi ATGC bên dưới
                    continue
                
                # Ghi tiêu đề vào file FASTA mới (EMU yêu cầu tiêu đề Fasta trùng khớp cột 1 của TSV)
                f_fasta.write(f">{current_acc}\n")
                
                # Ghi vào file taxonomy.tsv nếu TaxID này chưa từng được xử lý
                if current_acc not in processed_taxids:
                    try:
                        lineage_ids = ncbi.get_lineage(tax_id_int)
                        ranks = ncbi.get_rank(lineage_ids)
                        names = ncbi.get_taxid_translator(lineage_ids)
                        
                        target_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
                        lineage_dict = {rank: "" for rank in target_ranks}
                        
                        for t_id, rank in ranks.items():
                            if rank in target_ranks:
                                lineage_dict[rank] = names[t_id]
                        
                        row = [
                            current_acc,
                            lineage_dict['species'],
                            lineage_dict['genus'],
                            lineage_dict['family'],
                            lineage_dict['order'],
                            lineage_dict['class'],
                            lineage_dict['phylum'],
                            lineage_dict['superkingdom']
                        ]
                        f_tsv.write("\t".join(row) + "\n")
                        processed_taxids.add(current_acc)
                        valid_count += 1
                    except Exception:
                        pass
            else:
                # Ghi trình tự nucleotide nếu dòng header hợp lệ
                if current_acc:
                    f_fasta.write(f"{line}\n")

    print(f"\n🎉 Hoàn thành rực rỡ!")
    print(f"✅ Đã tạo bộ Database chuẩn EMU gồm {valid_count} loài phân loại sạch sẽ.")

if __name__ == "__main__":
    INPUT_FASTA = "raw_sequences.fasta"
    OUTPUT_FASTA = "species_taxid.fasta"
    OUTPUT_TSV = "taxonomy.tsv"
    
    if not os.path.exists(INPUT_FASTA):
        print(f"Không tìm thấy file {INPUT_FASTA} trong thư mục hiện tại!")
    else:
        build_emu_database_files = build_emu_database_smart(INPUT_FASTA, OUTPUT_FASTA, OUTPUT_TSV)