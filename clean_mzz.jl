# a script to manage files


# paths
root_dir = "/Users/ivancasas/GitHub/MIRS/Data/phd_data"
backup_dir = "/Users/ivancasas/GitHub/MIRS/Data/backup"
mkpath(backup_dir)


mzz_count = 0
opu_count = 0
txt_count = 0
otr_count = 0

mzz_list = []
opu_list = []
txt_list = []
otr_list = []


# count files, log them
for (root, _, files) in walkdir(root_dir)
    for file in files
        if endswith(file, ".mzz")
            mzz_count += 1
            push!(mzz_list, joinpath(root, file))

        elseif occursin(".txt", file)
            txt_count += 1
            push!(txt_list, joinpath(root, file))

        elseif occursin(r"\.\d+$", file) # any number of digits
            opu_count += 1
            push!(opu_list, joinpath(root, file))

        else
            otr_count += 1
            push!(otr_list, joinpath(root, file))
        end
    end
end

mzz_count
opu_count
txt_count
otr_count

mzz_list
opu_list
txt_list
otr_list

# create backup - rerun with care!
# for (root, _, files) in walkdir(root_dir)
#     for file in files
#         if occursin(r"\.\d+$", file) # any number of digits
#             cp(joinpath(root, file), joinpath(backup_dir, file), force = true)
#         end
#     end
# end

# check:
# backup_check = 0
# for (root, _, files) in walkdir(backup_dir)
#     for file in files
#         if occursin(r"\.\d+$", file) 
#             backup_check += 1
#         end
#     end
# end
# backup_check

# delete mzz and txt
for (root, _, files) in walkdir(root_dir)
    for file in files
        if endswith(file, ".mzz")
            rm(joinpath(root, file))
        elseif endswith(file, ".txt")
            rm(joinpath(root, file))
        end
    end
end

# deal with others

for i in otr_list
    println(i)
    println(isfile(i))
end

# common annoying files: ._, DS_Store
for (root, _, files) in walkdir(root_dir)
    for file in files
        if startswith(file, "._")
            rm(joinpath(root, file))
        elseif startswith(file, ".DS_Store")
            rm(joinpath(root, file))
        end
    end
end