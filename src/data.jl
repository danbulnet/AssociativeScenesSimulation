export ANAKG_SAMPLE_MAT_FILE, anakgload

using MAT

const ANAKG_SAMPLE_MAT_FILE = "data/ANAKG_data4.mat"

function anakgload(file=ANAKG_SAMPLE_MAT_FILE)::Matrix{Float64}
    matread(file)["ANAKG_graph"]
end