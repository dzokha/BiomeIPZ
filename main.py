from fastapi import FastAPI, BackgroundTasks
from pydantic import BaseModel
import subprocess
import os
import sys
import traceback

app = FastAPI(title="BiomeIPZ - Microbiome Pipeline API")

class PipelineRequest(BaseModel):
    job_id: str
    input_folder: str   # Thư mục chứa 36 file .fastq và metadata.csv
    output_folder: str  # Thư mục lưu kết quả

def run_bioinformatics_pipeline(job_id: str, input_dir: str, output_dir: str):
    """Hàm chạy ngầm gọi script pipeline chính"""
    log_file = os.path.join(output_dir, f"{job_id}_pipeline.log")
    os.makedirs(output_dir, exist_ok=True)
    
    with open(log_file, "w") as log:
        try:
            log.write(f"Đang khởi chạy pipeline bằng Python: {sys.executable}\n")
            log.write(f"Thư mục input: {input_dir}\n")
            log.write(f"Thư mục output: {output_dir}\n")
            log.write("-" * 50 + "\n")
            log.flush() # Ép ghi log ngay lập tức
            
            # Sử dụng sys.executable thay vì chữ "python"
            subprocess.run(
                [sys.executable, "pipeline.py", input_dir, output_dir],
                stdout=log, stderr=subprocess.STDOUT, check=True
            )
        except subprocess.CalledProcessError as e:
            log.write(f"\n[LỖI TIẾN TRÌNH] Pipeline trả về mã lỗi: {e.returncode}\n")
        except Exception as e:
            log.write(f"\n[LỖI HỆ THỐNG] Không thể khởi chạy tiến trình:\n{str(e)}\n")
            log.write(traceback.format_exc())

@app.post("/api/v1/run-pipeline")
async def trigger_pipeline(request: PipelineRequest, background_tasks: BackgroundTasks):
    # Kiểm tra thư mục đầu vào có tồn tại không
    if not os.path.exists(request.input_folder):
        return {"status": "error", "message": "Input folder not found!"}
    
    # Giao việc cho Background Task chạy ngầm
    background_tasks.add_task(
        run_bioinformatics_pipeline, 
        request.job_id, 
        request.input_folder, 
        request.output_folder
    )
    
    return {
        "status": "success", 
        "message": f"Job {request.job_id} is running in background.",
        "log_path": f"{request.output_folder}/{request.job_id}_pipeline.log"
    }

# Khởi chạy server: uvicorn main:app --host 0.0.0.0 --port 8000