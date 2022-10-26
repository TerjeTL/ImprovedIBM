function mesh = load_steady_state_solution(file, mesh_level_str, replace_0)
    %dir = "/solutions/" + mesh_level_str + "/time_dict";
    %time_list = h5read("../export_data.h5", dir);
    %str_time_list = compose('%0.6f', time_list);
    dir = "/solutions/" + mesh_level_str + "/steady_state/solution";
    mesh = h5read(file, dir);
    mesh(mesh==0) = replace_0;
end