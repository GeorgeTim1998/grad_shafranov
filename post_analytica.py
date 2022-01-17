import funcs as fu
from imports import *

print(colored("Date_Time is: %s" % fu.Time_name(), 'cyan'))
# vs_def_mesh_name = 'point_source'
# vs_square_size_name = vs_def_mesh_name
# # fu.Plot_umax_vs_def_mesh(vs_def_mesh_name)
# # fu.Plot_umax_vs_square_size(vs_square_size_name, default_mesh_size=fu.DEFAULT_MESH)

# fu.Plot_2D_data_together()
# это использовалось, чтобы строить графики сечения решения ГШ в случае точечного источника

fu.D_config(100)