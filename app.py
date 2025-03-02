from flask import Flask, request, jsonify
from flask_cors import CORS
import requests
import os
from dotenv import load_dotenv
import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Import the function from FinalMatrixBuilding
from FinalMatrixBuilding import matrix  # Adjust the function name as needed

app = Flask(__name__)
CORS(app)  # Enable communication between React frontend and Flask backend

# Load environment variables from .env file
load_dotenv()
deepseek_api_key = os.getenv("DEEPSEEK_API_KEY")

if not deepseek_api_key:
    raise ValueError("DEEPSEEK_API_KEY not found! Ensure it's correctly set in the .env file.")

# Directory where uploaded genome files will be stored
UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

# Load the gene-disease matrix from a CSV file (if available)
MATRIX_FILE = "cleaned_gene_disease_protein_mapping_final.csv"
matrix_data = None

if os.path.exists(MATRIX_FILE):
    matrix_data = pd.read_csv(MATRIX_FILE)

@app.route("/api/disease", methods=["POST"])
def get_disease_description():
    data = request.json
    disease = data.get("disease", "").strip()

    if not disease:
        return jsonify({"error": "Disease name required"}), 400

    try:
        # Prepare the request payload for DeepSeek API
        payload = {
            "model": "deepseek-chat",  # Replace with the appropriate DeepSeek model
            "messages": [
                {"role": "system", "content": "You are an expert in medical diseases."},
                {"role": "user", "content": f"Give a concise 50-word max description of the disease: {disease}"},
            ],
            "max_tokens": 100,
        }

        # Send request to DeepSeek API
        headers = {
            "Authorization": f"Bearer {deepseek_api_key}",
            "Content-Type": "application/json",
        }
        response = requests.post(
            "https://api.deepseek.com/v1/chat/completions",  # Replace with the actual DeepSeek API endpoint
            headers=headers,
            json=payload,
        )

        # Log full API response
        print("DeepSeek Response:", response.json())

        # Check for errors in the response
        if response.status_code != 200:
            return jsonify({"error": f"DeepSeek API error: {response.text}"}), response.status_code

        # Extract description
        response_data = response.json()
        if response_data.get("choices") and len(response_data["choices"]) > 0:
            description = response_data["choices"][0]["message"]["content"].strip()
            return jsonify({"description": description})

        return jsonify({"error": "No description found in API response"}), 500

    except requests.exceptions.RequestException as e:
        print("DeepSeek API Request Error:", str(e))
        return jsonify({"error": "Failed to connect to DeepSeek API"}), 500
    except Exception as e:
        print("General Exception:", str(e))
        return jsonify({"error": f"Internal Server Error: {str(e)}"}), 500

@app.route("/")
def home():
    return "Flask server is running!"

@app.route("/upload-genome", methods=["POST"])
def upload_genome():
    if "file" not in request.files:
        return jsonify({"message": "No file uploaded!"}), 400

    file = request.files["file"]

    # Ensure the file has a valid FASTA format
    if not file.filename.endswith((".fasta", ".fa", ".fna")):
        return jsonify({"message": "Invalid file type. Upload a FASTA genome file."}), 400

    # Save file
    file_path = os.path.join(UPLOAD_FOLDER, file.filename)
    file.save(file_path)

    return jsonify({"message": f"Genome file '{file.filename}' uploaded successfully!"})

@app.route("/upload-gff", methods=["POST"])
def upload_gff():
    if "file" not in request.files:
        return jsonify({"message": "No file uploaded!"}), 400

    file = request.files["file"]

    # Ensure the file has a valid GFF format
    if not file.filename.endswith(".gff"):
        return jsonify({"message": "Invalid file type. Upload a GFF file."}), 400

    # Save file
    file_path = os.path.join(UPLOAD_FOLDER, file.filename)
    file.save(file_path)

    return jsonify({"message": f"GFF file '{file.filename}' uploaded successfully!"})

@app.route("/process-files", methods=["POST"])
def process_files():
    if "gff_file" not in request.files or "fasta_file" not in request.files:
        return jsonify({"error": "Both GFF and FASTA files are required!"}), 400

    gff_file = request.files["gff_file"]
    fasta_file = request.files["fasta_file"]

    # Save files
    gff_path = os.path.join(UPLOAD_FOLDER, gff_file.filename)
    fasta_path = os.path.join(UPLOAD_FOLDER, fasta_file.filename)
    gff_file.save(gff_path)
    fasta_file.save(fasta_path)

    # Get the disease from the form data
    disease = request.form.get("disease", "").strip()
    if not disease:
        return jsonify({"error": "Disease name is required!"}), 400

    # Call the build_final_matrix function
    try:
        matrix = matrix(disease, gff_path, fasta_path)
        return jsonify({"message": "Files processed successfully!", "matrix": matrix})
    except Exception as e:
        print("Error processing files:", str(e))
        return jsonify({"error": f"Error processing files: {str(e)}"}), 500

@app.route("/api/matrix", methods=["GET"])
def get_matrix():
    if matrix_data is None:
        return jsonify({"error": "Matrix data not found. Ensure the CSV file is available."}), 404

    matrix_json = matrix_data.to_dict(orient="records")
    print("Matrix data fetched successfully.")
    return jsonify(matrix_json)

if __name__ == "__main__":
    app.run(debug=True)