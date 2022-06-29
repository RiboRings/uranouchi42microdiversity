using CSV, DataFrames

sample_files = readdir(joinpath(@__DIR__, "data"))
sample_files = sample_files[occursin.("microdiversity_summary", sample_files)]

genome_list = CSV.File("tmp.txt") |> DataFrame
genome_list[!, :selgen] = map(x -> string(x[1], "_matabat2bin.", x[2], ".fa"), split.(genome_list.selgen, "_"))

for cur_sample in sample_files

    cur_df = CSV.File(cur_sample) |> DataFrame

    keep_rows = map(x -> x ∈ genome_list.selgen, cur_df.genome)
    out_df = cur_df[keep_rows, :]
    
    CSV.write(string(@__DIR__, "/data/", cur_sample, "_short.csv"), out_df)

end
