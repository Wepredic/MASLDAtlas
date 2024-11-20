
# MASLDAtlas: Single-cell Data Analysis and Visualization

MASLDAtlas is a user-friendly and interactive R Shiny application designed for the analysis and visualization of single-cell data. 



## Installation

### Prerequisites

Before running the application, ensure the following are installed:

1. **R (≥ 4.0.0)**: [Download R](https://cran.r-project.org/)
2. **RStudio (Optional)**: [Download RStudio](https://posit.co/download/rstudio-desktop/)
3. **Docker & Docker Compose**: Install Docker and Docker Compose for containerized deployment.

### Clone the Repository

```bash
git clone https://github.com/Wepredic/MASLDAtlas.git
cd MASLDAtlas
```

---

## Running the Application

### Option 1: Run Locally in R

1. Install required R packages:
    ```R
    Rscript ./data/installpackages.R
    Rscript ./data/reticulate_create_env.R
    ```
2. Run the following command to launch the app:
    ```R
    Rscript ./data/app.R
    ```

The application will launch in your default web browser.

---

### Option 2: Run via Docker

#### Build the Docker Image

To build the Docker image, run:

```bash
docker build -t masldatlas .
```

#### Deploy with Docker Compose

Using the provided `docker-compose.yml`, you can easily set up the application along with its dependencies. To start the application, run:

```bash
docker-compose up
```

The app will be accessible at `http://localhost:6868`.

---

## Docker Architecture

### Dockerfile
The `Dockerfile` defines the application environment:
- **Base Image**: R-based image for Shiny applications.
- **Dependencies**: Installs all necessary R packages.
- **Application**: Adds the R Shiny app to the container.

### Docker Compose
The `docker-compose.yml` includes:
- **Volumes**: Maps the application directory to the container for seamless updates.
- **Ports**: Exposes the application on port `6868`.
- **Services**:
  - `app`: The R Shiny server hosting the MASLDAtlas application.

---

## Data Input

- **Accepted Formats**: Loom files are the primary format for single-cell data.
- **Usage**: Upload Loom files via the app interface or place them in the designated `data` directory (when using Docker).

---

## Contributions

We welcome contributions! Please submit an issue or pull request for bugs, feature requests, or suggestions.


