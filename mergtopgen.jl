using CSV, DataFrames

sample_files = readdir(joinpath(@__DIR__, "data"))
sample_files = sample_files[occursin.("short", sample_files)]

raw_sample_names = replace.(sample_files, "_snvs_summary.csv_short.csv" => "")
sample_names = map(x -> string(x[4:end], x[1:3]), raw_sample_names)

println("Importing files")

sample_file_list = map(x -> CSV.File(joinpath(@__DIR__, "data", x)) |> DataFrame, sample_files)

println("Concatenating files")

snvs_df = sample_file_list[1]
snvs_df[!, :sample] .= sample_names[1]

for i in 2:length(sample_file_list)

    cur_df = sample_file_list[i]
    cur_df[!, :sample] .= sample_names[i]

    global snvs_df = vcat(snvs_df, cur_df)

end

println("Writing output")

CSV.write(joinpath(@__DIR__, "data/snvs.csv"), snvs_df)
