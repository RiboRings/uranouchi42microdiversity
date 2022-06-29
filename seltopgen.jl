using CSV, DataFrames

sample_files = readdir(joinpath(@__DIR__, "data"))
sample_files = sample_files[occursin.("snvs_summary", sample_files)]

genome_list = CSV.File("top_genome_list.txt") |> DataFrame
genome_list[!, :selgen] = map(x -> string(x[1], "_matabat2bin.", x[2], ".fa"), split.(genome_list.selgen, "_"))

for cur_sample in sample_files

    println(cur_sample)
    cur_df = CSV.File(joinpath(@__DIR__, "data", cur_sample)) |> DataFrame

    keep_rows = map(x -> x ∈ genome_list.selgen, cur_df.genome)
    out_df = cur_df[keep_rows, :]
    
    CSV.write(string(@__DIR__, "/data/", cur_sample, "_short.csv"), out_df)
    println(string(cur_sample, " done"))

end
