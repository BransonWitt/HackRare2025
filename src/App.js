import React, { useState, useEffect } from "react";
import "./App.css";
import logo from "./GenieBottleOfficial.png"; 
import axios from "axios";

const App = () => {
  const [disease, setDisease] = useState("");
  const [description, setDescription] = useState("");
  const [matrix, setMatrix] = useState([]);
  const [genes, setGenes] = useState([]);
  const [file, setFile] = useState(null);
  const [filePreview, setFilePreview] = useState("");
  const [message, setMessage] = useState("");
  const [diseasesList, setDiseasesList] = useState([]); // Stores diseases from CSV
  const [filteredDiseases, setFilteredDiseases] = useState([]); // Filtered list for dropdown
  const [gffFile, setGffFile] = useState(null);
  const [gffFilePreview, setGffFilePreview] = useState("");
  const [gffMessage, setGffMessage] = useState("");

  useEffect(() => {
    fetch("/cleaned_gene_disease_protein_mapping_final.csv")
      .then((response) => response.text())
      .then((data) => {
        const rows = data.split("\n").map((row) => row.trim().split(",")); // Split into rows & columns
        const diseaseLabels = rows.slice(1).map((row) => row[1]); // Extract Disease Label column
        const uniqueDiseases = [...new Set(diseaseLabels.filter(Boolean))];

        setDiseasesList(uniqueDiseases);
        setFilteredDiseases(uniqueDiseases);
      })
      .catch((error) => console.error("Error loading diseases:", error));
  }, []);

  // Filter diseases when typing in the search bar
  const handleSearchChange = (event) => {
    const query = event.target.value.toLowerCase();
    setDisease(query);

    const filtered = diseasesList.filter((d) =>
      d.toLowerCase().includes(query)
    );
    setFilteredDiseases(filtered);
  };

  // Select a disease from the dropdown
  const handleSelectDisease = (selectedDisease) => {
    setDisease(selectedDisease);
    setFilteredDiseases([]);
  };

  // handle disease search
  const handleSearch = async () => {
    if (!disease) return alert("Please select a disease.");
  
    try {
      const response = await fetch("http://127.0.0.1:5000/api/disease", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({ disease }),
      });
  
      const data = await response.json();
      if (data.description) {
        setDescription(data.description);
      } else {
        setDescription("No description found.");
      }
    } catch (error) {
      console.error("Error fetching disease description:", error);
      setDescription("Error fetching disease information.");
    }
  };

  // Function to handle genome file upload
  const handleFileUpload = (event) => {
    const uploadedFile = event.target.files[0];

    if (uploadedFile) {
      if (uploadedFile.name.endsWith(".fasta") || uploadedFile.name.endsWith(".fa") || uploadedFile.name.endsWith(".fna")) {
        setFile(uploadedFile);

        const reader = new FileReader();
        reader.onload = (e) => {
          const fileContent = e.target.result.trim();
          
          if (fileContent.startsWith(">")) {
            setMessage("Valid FASTA genome file uploaded.");
          } else {
            setMessage("Invalid FASTA file: Missing '>' header.");
            setFile(null);
          }
        };
        reader.readAsText(uploadedFile);
      } else {
        setMessage("Please upload a valid FASTA file (.fasta or .fa).");
        setFile(null);
      }
    }
  };

  // send the genome file to Flask
  const handleSubmit = async () => {
    if (!file) {
      alert("Please upload a valid genome FASTA file.");
      return;
    }

    const formData = new FormData();
    formData.append("file", file);

    try {
      const response = await fetch("http://127.0.0.1:5000/upload-genome", {
        method: "POST",
        body: formData,
      });

      const result = await response.json();
      setMessage(result.message);
    } catch (error) {
      console.error("Error uploading genome:", error);
      setMessage("Error uploading genome file.");
    }
  };

  // Function to handle GFF file upload
  const handleGffFileUpload = (event) => {
    const uploadedFile = event.target.files[0];

    if (uploadedFile) {
      if (uploadedFile.name.endsWith(".gff")) {
        setGffFile(uploadedFile);

        const reader = new FileReader();
        reader.onload = (e) => {
          const fileContent = e.target.result.trim();
          setGffMessage("GFF file uploaded.");
        };
        reader.readAsText(uploadedFile);
      } else {
        setGffMessage("Please upload a valid GFF file (.gff).");
        setGffFile(null);
      }
    }
  };

  // send the GFF file to Flask
  const handleGffSubmit = async () => {
    if (!gffFile) {
      alert("Please upload a valid GFF file.");
      return;
    }

    const formData = new FormData();
    formData.append("file", gffFile);

    try {
      const response = await fetch("http://127.0.0.1:5000/upload-gff", {
        method: "POST",
        body: formData,
      });

      const result = await response.json();
      setGffMessage(result.message);
    } catch (error) {
      console.error("Error uploading GFF file:", error);
      setGffMessage("Error uploading GFF file.");
    }
  };
  const handleProcessFiles = async () => {
    if (!gffFile || !file) {
      alert("Please upload both GFF and FASTA files.");
      return;
    }
  
    const formData = new FormData();
    formData.append("gff_file", gffFile);
    formData.append("fasta_file", file);
  
    try {
      const response = await fetch("http://127.0.0.1:5000/process-files", {
        method: "POST",
        body: formData,
      });
  
      if (result.matrix && Array.isArray(result.matrix)) {
        setMatrix(result.matrix);  // Ensure matrix is stored in state
      } else {
        console.error("Invalid matrix format received.");
        setMatrix([]); // Reset matrix if response is invalid
      }
    } catch (error) {
      console.error("Error processing files:", error);
      setMatrix([]); // Reset on error
    }
  };

  return (
    <div>
      {/* Navigation Bar */}
      <nav className="navbar">
        <div className="logo-container">
          <img src={logo} alt="Gene-ie Logo" className="logo" />
          <span className="animated-text">Gene-ie</span> {/* Animated text */}
        </div>
        <div className="nav-links">
          <button className="nav-button">Home</button>
          <button className="nav-button">About</button>
          <button className="nav-button">Contact</button>
        </div>
      </nav>

      {/* Hero Section */}
      <header className="hero">
        <h1>Gene-ie</h1>
        <p>Discover genetic insights and potential treatments with our powerful search tool.</p>

        {/* Searchable Dropdown */}
        <div className="search-section">
          <input
            type="text"
            placeholder="Search for a disease..."
            value={disease}
            onChange={handleSearchChange}
            className="search-input"
          />
          <button onClick={handleSearch}>Search</button>

          {/* Dropdown List */}
          {disease && filteredDiseases.length > 0 && (
            <ul className="dropdown">
              {filteredDiseases.slice(0, 5).map((d, idx) => (
                <li key={idx} onClick={() => handleSelectDisease(d)}>
                  {d}
                </li>
              ))}
            </ul>
          )}
        </div>
      </header>

      {/* Upload Section for FASTA files */}
      <div className="upload-section">
        <h2>Upload Your Genome (FASTA Format)</h2>
        <input type="file" accept=".fna,.fasta,.fa" onChange={handleFileUpload} />
        {file && <p>Uploaded: {file.name}</p>}
        {message && <p className="upload-message">{message}</p>}
        {filePreview && (
          <div className="file-preview">
            <pre>{filePreview}</pre>
          </div>
        )}
        <button onClick={handleSubmit} className="upload-button" disabled={!file}>
        Submit Genome
        </button>
    </div>

      {/* Upload Section for GFF files */}
      <div className="upload-section">
        <h2>Upload Your GFF File</h2>
        <input type="file" accept=".gff" onChange={handleGffFileUpload} />
        {gffFile && <p>Uploaded: {gffFile.name}</p>}
        {gffMessage && <p className="upload-message">{gffMessage}</p>}
        {gffFilePreview && (
          <div className="file-preview">
            <pre>{gffFilePreview}</pre>
          </div>
        )}
        <button onClick={handleGffSubmit} className="upload-button" disabled={!gffFile}>
          Submit GFF File
        </button>
      </div>

      {/* Description & Matrix */}
      <div className="info-container">
        <div className="description">
          <h2>Description</h2>
          <p>{description}</p>
        </div>

        <div className="matrix">
          <h2>Gene-Disease Matrix</h2>
          {matrix.length > 0 ? (
            <table>
              <thead>
                <tr>
                  {matrix[0].map((col, idx) => (
                    <th key={idx}>{col}</th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {matrix.slice(1).map((row, idx) => (
                  <tr key={idx}>
                    {row.map((cell, i) => (
                      <td key={i}>{cell}</td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          ) : (
            <p>exampleMatrix.csv</p>
          )}
        </div>
      </div>

      {/* Buttons for Pathway & Medicine */}
      <div className="buttons">
        <a
          href={`https://reactome.org/content/query?q=${disease}`}
          target="_blank"
          rel="noopener noreferrer"
        >
          <button className="pathway-button">Pathway</button>
        </a>
        <a
          href={`https://go.drugbank.com/unearth/q?search=${disease}&searcher=drugs`}
          target="_blank"
          rel="noopener noreferrer"
        >
          <button className="medicine-button">Medicine/Drugs</button>
        </a>
      </div>

      {/* Footer */}
      <footer className="footer">
        <p>© 2025 Gene-ie | Advancing Genomic Research</p>
      </footer>
    </div>
  );
};
export default App;
