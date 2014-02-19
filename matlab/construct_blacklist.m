function blacklist = construct_blacklist(whitelist,full_list)

blacklist = full_list(~ismember(full_list,whitelist));