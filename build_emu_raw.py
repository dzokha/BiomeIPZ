import os
import sys
import gzip
from ete3 import NCBITaxa

def build_emu_database_from_catalog(fasta_in, catalog_gz, fasta_out, tsv_out):
    print("🌳 Bước 1: Khởi tạo dữ liệu cây phả hệ sinh giới...")
    ncbi = NCBITaxa()
    
    print("🔍 Bước 2: Quét file FASTA để lấy danh sách mã Accession cần tìm...")
    target_accs = set()
    with open(fasta_in, 'r') as f:
        for line in f:
            if line.startswith(">"):
                # Lấy mã Accession chuẩn, ví dụ: "NR_113675.1"
                acc = line[1:].split()[0]
                target_accs.add(acc)
                
    total_targets = len(target_accs)
    print(f"✅ Tìm thấy {total_targets} mã Accession cần ánh xạ.")

    print("\n🗺️ Bước 3: Nạp danh bạ RefSeq Catalog và ánh xạ ngoại tuyến...")
    acc_to_taxid = {}
    
    # Đọc trực tiếp file catalog nén (.gz) theo từng dòng để tiết kiệm RAM
    with gzip.open(catalog_gz, 'rt') as f:
        for line in f:
            # File catalog phân tách bằng dấu tab
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                taxid_str = parts[0]  # Cột 1 là TaxID dạng số
                acc_ver = parts[3]    # Cột 4 là Accession.version (e.g., NR_113675.1)
                
                if acc_ver in target_accs:
                    acc_to_taxid[acc_ver] = int(taxid_str)
                    
                    if len(acc_to_taxid) % 500 == 0:
                        sys.stdout.write(f"\r⏳ Đã khớp mã thành công: {len(acc_to_taxid)}/{total_targets} chủng loài...")
                        sys.stdout.flush()

    print(f"\n✅ Ánh xạ thành công! Đã tìm thấy TaxID cho {len(acc_to_taxid)} chủng loài.")

    print("\n🧬 Bước 4: Đúc file species_taxid.fasta và taxonomy.tsv đúng chuẩn EMU...")
    with open(fasta_in, 'r') as f_in, \
         open(fasta_out, 'w') as f_fasta, \
         open(tsv_out, 'w') as f_tsv:
        
        # Ghi Header chuẩn của EMU
        f_tsv.write("tax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tclade\n")
        
        processed_accs = set()
        current_acc = None
        valid_count = 0

        for line in f_in:
            if line.startswith(">"):
                current_acc = line[1:].split()[0]
                
                tax_id_int = acc_to_taxid.get(current_acc)
                if not tax_id_int:
                    current_acc = None  # Bỏ qua nếu chủng này không có TaxID trong danh bạ
                    continue
                
                f_fasta.write(f">{current_acc}\n")
                
                if current_acc not in processed_accs:
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
                        processed_accs.add(current_acc)
                        valid_count += 1
                    except Exception:
                        pass
            else:
                if current_acc:
                    f_fasta.write(f"{line}\n")

    print(f"\n🎉 HOÀN THÀNH XUẤT SẮC!")
    print(f"✅ Đã tạo bộ dữ liệu sạch cho EMU gồm {valid_count} loài phân loại hoàn chỉnh.")

if __name__ == "__main__":
    INPUT_FASTA = "sequence.fasta"   # File 41.9 MB tải từ NCBI
    
    # Tự động tìm file Catalog trong thư mục (bạn đỡ phải đổi tên file khi tải bản mới)
    catalog_files = [f for f in os.listdir('.') if f.startswith('RefSeq-release') and f.endswith('.catalog.gz')]
    
    if not os.path.exists(INPUT_FASTA):
        print(f"❌ Không tìm thấy file {INPUT_FASTA}. Hãy để file 41.9MB vào đây.")
    elif not catalog_files:
        print(f"❌ Không tìm thấy file danh bạ dạng 'RefSeq-releaseXXXX.catalog.gz'. Hãy tải nó về và để cạnh file script này.")
    else:
        CATALOG_GZ = catalog_files[0]
        print(f"📦 Tìm thấy file danh bạ: {CATALOG_GZ}")
        OUTPUT_FASTA = "species_taxid.fasta"
        OUTPUT_TSV = "taxonomy.tsv"
        build_emu_database_from_catalog(INPUT_FASTA, CATALOG_GZ, OUTPUT_FASTA, OUTPUT_TSV)