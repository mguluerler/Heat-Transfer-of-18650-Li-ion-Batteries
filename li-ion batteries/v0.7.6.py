# Version Notes(v0.7.6):    1)For finding of area in dT/dz, pi*r^2 added!
#                           2)T_mean started calculating with area!
#                           3)k_z and k_r are not same values now!
#                           4)Convection_type added for nodeFormulas class, other convection types could be used after formulations added
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
#import numba
#import translater_outfile

T_initial, max_flowtime, t_part, z_part, r_part, dz, dr, dt, t, T, g, density_Battery, specHeat_Battery, alfa, \
h_nearused, h_topused, h_bottomused, img = [None]*18
class solvingParams():
    global T_initial, max_flowtime, t_part, z_part, r_part, dz, dr, dt, t, T, g
    T_initial = 22 + 273.15
    max_flowtime = 300 #s
    t_part, z_part, r_part = 30000, 65, 9
    dz = 0.065/z_part #m
    dr = 0.009/r_part #m
    dt = max_flowtime/t_part
    T = np.zeros((t_part+1,z_part+1, r_part+1), dtype=float) + (T_initial+0.01)
    t = 0
    g = 9.81
class BattProp():
    global density_Battery, specHeat_Battery, alfa
    density_Battery = 1694.071
    specHeat_Battery = 903.4735
    alfa = 1 /(density_Battery*specHeat_Battery)
class dataSavers():
    global h_nearused, h_topused, h_bottomused
    h_nearused = []
    h_topused = []
    h_bottomused = []
    for emptylistopener in range(0, math.floor(t_part)):
        h_nearused.append([])
        h_topused.append([])
        h_bottomused.append([])

# In one dimension
#         T, density, specific heat, thermal conductivity, thermal diffusivity, dynamic viscosity, kinematic viscosity, prandtl number
thermotable_air = [[ 0,   1.292,          1006,              0.02364,    1.818*(10**(-5)),  1.729*(10**(-5)),    1.338*(10**(-5)),        0.7362],[
          5,   1.269,          1006,              0.02401,    1.880*(10**(-5)),  1.754*(10**(-5)),    1.382*(10**(-5)),        0.7350],[
         10,   1.246,          1006,              0.02439,    1.944*(10**(-5)),  1.778*(10**(-5)),    1.426*(10**(-5)),        0.7336],[
         15,   1.225,          1007,              0.02476,    2.009*(10**(-5)),  1.802*(10**(-5)),    1.470*(10**(-5)),        0.7323],[
         20,   1.204,          1007,              0.02514,    2.074*(10**(-5)),  1.825*(10**(-5)),    1.516*(10**(-5)),        0.7309],[
         25,   1.184,          1007,              0.02551,    2.141*(10**(-5)),  1.849*(10**(-5)),    1.562*(10**(-5)),        0.7296],[
         30,   1.164,          1007,              0.02588,    2.208*(10**(-5)),  1.872*(10**(-5)),    1.608*(10**(-5)),        0.7282],[
         35,   1.145,          1007,              0.02625,    2.277*(10**(-5)),  1.895*(10**(-5)),    1.655*(10**(-5)),        0.7268],[
         40,   1.127,          1007,              0.02662,    2.346*(10**(-5)),  1.918*(10**(-5)),    1.702*(10**(-5)),        0.7255],[
         45,   1.109,          1007,              0.02699,    2.416*(10**(-5)),  1.941*(10**(-5)),    1.750*(10**(-5)),        0.7241],[
         50,   1.092,          1007,              0.02735,    2.487*(10**(-5)),  1.963*(10**(-5)),    1.798*(10**(-5)),        0.7228],[
         60,   1.059,          1007,              0.02808,    2.632*(10**(-5)),  2.008*(10**(-5)),    1.896*(10**(-5)),        0.7202],[
         70,   1.028,          1007,              0.02881,    2.780*(10**(-5)),  2.052*(10**(-5)),    1.995*(10**(-5)),        0.7177]]
thermotabletempscelc_air =  [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0]
thermotabletempskelvin_air = list(map(lambda kelvinadder: kelvinadder + 273.15, thermotabletempscelc_air))

class materials:
        def __init__(self, T_):
                self.T = T_
                self.Tf = self.T
                self.B = 1/self.Tf
                upper_index = thermotabletempskelvin_air.index(min(list(filter(lambda sicak: sicak > self.T, thermotabletempskelvin_air))))
                lower_index = thermotabletempskelvin_air.index(max(list(filter(lambda sicak: sicak < self.T, thermotabletempskelvin_air))))
                oran = (self.Tf - thermotable_air[lower_index][0]) / (thermotable_air[upper_index][0] - thermotable_air[lower_index][0])
                self.k = (thermotable_air[upper_index][3] - thermotable_air[lower_index][3]) * oran + thermotable_air[lower_index][3]
                self.Pr = (thermotable_air[upper_index][7] - thermotable_air[lower_index][7]) * oran + thermotable_air[lower_index][7]
                self.ro = (thermotable_air[upper_index][1] - thermotable_air[lower_index][1]) * oran + thermotable_air[lower_index][1]
                self.c = (thermotable_air[upper_index][2] - thermotable_air[lower_index][2]) * oran + thermotable_air[lower_index][2]
                self.kinV = (thermotable_air[upper_index][6] - thermotable_air[lower_index][6]) * oran + thermotable_air[lower_index][6]
        def CalculateNewCond(self):
                upper_index = thermotabletempskelvin_air.index(min(list(filter(lambda sicak: sicak > self.T, thermotabletempskelvin_air))))
                lower_index = thermotabletempskelvin_air.index(max(list(filter(lambda sicak: sicak < self.T, thermotabletempskelvin_air))))
                oran = (self.Tf - thermotable_air[lower_index][0]) /(thermotable_air[upper_index][0] - thermotable_air[lower_index][0])
                self.k = (thermotable_air[upper_index][3] - thermotable_air[lower_index][3])*oran + thermotable_air[lower_index][3]
                self.Pr = (thermotable_air[upper_index][7] - thermotable_air[lower_index][7])*oran + thermotable_air[lower_index][7]
                self.ro = (thermotable_air[upper_index][1] - thermotable_air[lower_index][1])*oran + thermotable_air[lower_index][1]
                self.c = (thermotable_air[upper_index][2] - thermotable_air[lower_index][2])*oran + thermotable_air[lower_index][2]
                self.kinV = (thermotable_air[upper_index][6] - thermotable_air[lower_index][6])*oran + thermotable_air[lower_index][6]
class k_find:
    k_r = 1.4685346
    k_z = 29.5753143
class h_find:
    @staticmethod
    def calculate(type):
        T_o = T[t][j][i]
        h_find.T_f(T_o)
        material_air.CalculateNewCond()
        if type == "near_natural":
            Gra_L = h_find.Gr_L(T_o, L=0.065)
            Ra = Gra_L*material_air.Pr
            nusselt = h_find.Nu_L_CylSide(Ra, L=0.065, d=0.018)
            h_finded = h_find.h(nusselt, L=0.065)
            h_nearused[t].append(h_finded)
            return h_finded
        elif type == "top_natural":
            Gra_L = h_find.Gr_L(T_o, L=0.018/4)
            Ra = Gra_L*material_air.Pr
            if Ra <= 8000000:
                nusselt = 0.54 * (Ra ** 0.25)
            else:
                nusselt = 0.15 * (Ra ** (1.0 / 3.0))
            h_finded = h_find.h(nusselt, L=0.018/4)
            h_topused[t].append(h_finded)
            return h_finded
        elif type == "bottom_natural":
            Gra_L = h_find.Gr_L(T_o, L=0.018/4)
            Ra = Gra_L*material_air.Pr
            nusselt = 0.27 * (Ra ** 0.25)
            h_finded = h_find.h(nusselt, L=0.018/4)
            h_bottomused[t].append(h_finded)
            return h_finded
    @staticmethod
    def T_f(T_o):
            material_air.Tf = (T_o + material_air.T)/2
            material_air.B = 1/material_air.Tf
    @staticmethod
    def Gr_L(T_o, L):
            return (g*material_air.B* (T_o-material_air.T) * (L**3)) / (material_air.kinV**2)
    @staticmethod
    def Nu_L_CylSide(Ra, L, d):
            return ((4/3) * (((7*material_air.Pr/5)/(20+(21*material_air.Pr)))**0.25) * (Ra**0.25)) + ((4/35) * (272+315*material_air.Pr) * L)/((64 + 63*material_air.Pr) * d)
    @staticmethod
    def h(Nu, L):
            return Nu*material_air.k/L
class R_find:
    @staticmethod
    def Area_finder(type, i):
        if type == "conv_near":
            Area = 2 * math.pi * i * dr * dz
            if j == 0 or j == z_part:
                return Area/2
            else:
                return Area

        elif type == "conv_top":
            Area = math.pi*(((i+1)*dr)**2 - ((i-1)*dr)**2)
            if i == r_part:
                return math.pi*(((i)*dr)**2 - ((i-1)*dr)**2)
            elif i == 0:
                Area = math.pi * ((dr / 2)**2)
                return Area
            else:
                return Area

        elif type == "cond_r":
            Area = 2 * math.pi * i * dr * dz
            if i == 0:
                Area = 2 * math.pi * (dr / 4) * dz
                if j == 0 or j == z_part:
                    return Area/2
                else:
                    return Area
            else:
                if j == 0 or j == z_part:
                    return Area/2
                else:
                    return Area

        elif type == "cond_z":
            Area = math.pi*(((i+1)*dr)**2 - ((i-1)*dr)**2)
            if i == r_part:
                return math.pi*(((i)*dr)**2 - ((i-1)*dr)**2)
            elif i == 0:
                Area = math.pi * ((dr / 2)**2)
                return Area
            else:
                return Area

    @staticmethod
    def dt_dr(type, i):
        if type == "forward":
            if i == 0:
                return dr / k_find.k_r
            else:
                return math.log((i+1)/i, math.e) * i * dr / k_find.k_r
        elif type == "backward":
            if i == 1:
                return dr / k_find.k_r
            else:
                return math.log(i/(i-1), math.e) * i * dr / k_find.k_r
    @staticmethod
    def dt_dz():
        return dz / k_find.k_z
def volume_node():
    if i == r_part and (j == 0 or j == z_part):
        gen_volume = 2 * math.pi * i * dr * (dr / 2) * (dz / 2)
    elif i == 0 and (j == 0 or j == z_part):
        gen_volume = math.pi * ((dr / 2)**2) * (dz / 2)
    elif i == r_part:
        gen_volume = 2 * math.pi * i * dr * (dr / 2) * dz
    elif i == 0:
        gen_volume = math.pi * ((dr / 2)**2) * dz
    elif j == 0 or j == z_part:
        gen_volume = 2 * math.pi * i * dr * dr * (dz / 2)
    else:
        gen_volume = 2 * math.pi * i * dr * dr * dz
    return gen_volume
def e_generation():
    return (-0.000000000386585918245514000000*((t*dt)**6)
            + 0.000000632699976906709000000000*((t*dt)**5)
            - 0.000376441425752549000000000000*((t*dt)**4)
            + 0.094957183016276800000000000000*((t*dt)**3)
            - 8.147904982120960000000000000000*((t*dt)**2)
            - 170.423649147152000000000000000000*((t*dt)**1)
            + 145193.154369295000000000000000000000) * volume_node()
class nodeFormulas:
    """
        convr => i == r_part
        convz => z == 0 or z == z_part
    """
    convection_type = ""
    @staticmethod
    def T_condall(t, j, i): # T[time][z][r]
        if i==0:
            T[t + 1][j][i] = (((T[t][j][i + 1] - T[t][j][i]) * R_find.Area_finder(type="cond_r", i=i) / R_find.dt_dr(type="forward", i=i) +
                               (T[t][j - 1][i] - T[t][j][i]) * R_find.Area_finder(type="cond_z", i=i) / R_find.dt_dz() +
                               (T[t][j + 1][i] - T[t][j][i]) * R_find.Area_finder(type="cond_z", i=i) / R_find.dt_dz() +
                               e_generation()) * dt * alfa / volume_node() +
                              T[t][j][i])

        else:
            T[t+1][j][i] = (((T[t][j][i + 1] - T[t][j][i]) * R_find.Area_finder(type="cond_r", i=i) / R_find.dt_dr(type="forward", i=i) +
                             (T[t][j][i - 1] - T[t][j][i]) * R_find.Area_finder(type="cond_r", i=i) / R_find.dt_dr(type="backward", i=i) +
                             (T[t][j - 1][i] - T[t][j][i]) * R_find.Area_finder(type="cond_z", i=i) / R_find.dt_dz() +
                             (T[t][j + 1][i] - T[t][j][i]) * R_find.Area_finder(type="cond_z", i=i) / R_find.dt_dz() +
                             e_generation()) * dt*alfa / volume_node() +
                            T[t][j][i])
    @staticmethod
    def T_convr(t, j, i):
        h_finded = h_find.calculate(type=("near"+nodeFormulas.convection_type))
        T[t + 1][j][i] = ((((T[t][j - 1][i] - T[t][j][i]) * R_find.Area_finder(type="cond_z", i=i) / R_find.dt_dz() +
                            (T[t][j + 1][i] - T[t][j][i]) * R_find.Area_finder(type="cond_z", i=i) / R_find.dt_dz()) +
                            (T[t][j][i - 1] - T[t][j][i]) * R_find.Area_finder(type="cond_r", i=i) / R_find.dt_dr(type="backward", i=i) +
                            (h_finded * (material_air.T - T[t][j][i]) * R_find.Area_finder(type="conv_near", i=i)) +
                           e_generation()) * dt * alfa / volume_node() +
                          T[t][j][i])
    @staticmethod
    def T_convz(t, j, i):
        if i == 0:
            eq_tempr = ((T[t][j][i + 1] - T[t][j][i]) * R_find.Area_finder(type="cond_r", i=i) / R_find.dt_dr(type="forward", i=i))
        else:
            eq_tempr = ((T[t][j][i + 1] - T[t][j][i]) * R_find.Area_finder(type="cond_r", i=i) / R_find.dt_dr(type="forward", i=i) +
                        (T[t][j][i - 1] - T[t][j][i]) * R_find.Area_finder(type="cond_r", i=i) / R_find.dt_dr(type="backward", i=i))
        if j == 0:
            h_finded = h_find.calculate(type=("top"+nodeFormulas.convection_type))
            T[t + 1][j][i] =((eq_tempr +  # dr
                            (T[t][j + 1][i] - T[t][j][i]) * R_find.Area_finder(type="cond_z", i=i) / R_find.dt_dz() +
                            (h_finded * (material_air.T - T[t][j][i]) * R_find.Area_finder(type="conv_top", i=i)) +
                            e_generation()) * alfa * dt / volume_node() +
                            T[t][j][i])
        else:
            h_finded = h_find.calculate(type=("bottom"+nodeFormulas.convection_type))
            T[t + 1][j][i] = (eq_tempr +  # dr
                             ((T[t][j - 1][i] - T[t][j][i]) * R_find.Area_finder(type="cond_z", i=i) / R_find.dt_dz() +
                              (h_finded * (material_air.T - T[t][j][i]) * R_find.Area_finder(type="conv_top", i=i)) +
                              e_generation()) * alfa * dt / volume_node() +
                              T[t][j][i])
    @staticmethod
    def T_convcorner(t, j, i):
        h_r = h_find.calculate(type=("near"+nodeFormulas.convection_type))
        if j == 0:
            h_z = h_find.calculate(type=("top"+nodeFormulas.convection_type))
            T[t + 1][j][i] =(((T[t][j + 1][i] - T[t][j][i]) * R_find.Area_finder(type="cond_z", i=i) / R_find.dt_dz() +
                              (T[t][j][i - 1] - T[t][j][i]) * R_find.Area_finder(type="cond_r", i=i) / R_find.dt_dr(type="backward", i=i) +
                              (h_r * (material_air.T - T[t][j][i]) * R_find.Area_finder(type="conv_near", i=i)) +
                              (h_z * (material_air.T - T[t][j][i]) * R_find.Area_finder(type="conv_top", i=i)) +
                              e_generation()) * dt*alfa / volume_node() +
                             T[t][j][i])
        else:
            h_z = h_find.calculate(type=("bottom"+nodeFormulas.convection_type))
            T[t + 1][j][i] = (((T[t][j - 1][i] - T[t][j][i]) * R_find.Area_finder(type="cond_z", i=i) / R_find.dt_dz() +
                               (T[t][j][i - 1] - T[t][j][i]) * R_find.Area_finder(type="cond_r", i=i) / R_find.dt_dr(type="backward", i=i) +
                               (h_r * (material_air.T - T[t][j][i]) * R_find.Area_finder(type="conv_near", i=i)) +
                               (h_z * (material_air.T - T[t][j][i]) * R_find.Area_finder(type="conv_top", i=i)) +
                               e_generation()) * dt * alfa / volume_node() +
                              T[t][j][i])
class graphs:
    scale = 10
    @staticmethod
    def tempmap_zr(time, ax):
        global img
        ax.set_title("Battery Temperature Map\n[Blue: Cold, Red: Hot]")
        scale = graphs.scale
        img = np.zeros(((z_part) * scale + 1, (r_part) * scale + 1, 3), np.uint8)
        xdraw = 0
        ydraw = 0
        while xdraw < r_part * scale:
            for x in range(xdraw, xdraw + scale):
                delta_a = x - xdraw
                if x == xdraw:
                    continue
                while ydraw < z_part * scale:
                    for y in range(ydraw, ydraw + scale):
                        delta_b = y - ydraw
                        if y == ydraw:
                            continue
                        base_color_deg = (T[time, math.floor(y / scale), math.floor(x / scale)] - np.min(T[time])) / (
                                np.max(T[time]) - np.min(T[time])) * 500
                        nexty_color_deg = (T[time, math.floor((y + scale) / scale), math.floor(x / scale)] - np.min(
                            T[time])) / (
                                                  np.max(T[time]) - np.min(T[time])) * 500
                        nextx_color_deg = (T[time, math.floor(y / scale), math.floor((x + scale) / scale)] - np.min(
                            T[time])) / (
                                                  np.max(T[time]) - np.min(T[time])) * 500
                        nextcorner_color_deg = (T[time, math.floor((y + scale) / scale), math.floor(
                            (x + scale) / scale)] - np.min(T[time])) / (
                                                       np.max(T[time]) - np.min(T[time])) * 500
                        a_color_deg = ((nextx_color_deg - base_color_deg) * delta_a / scale) + base_color_deg
                        b_color_deg = ((nextcorner_color_deg - nexty_color_deg) * delta_a / scale) + nexty_color_deg
                        color_deg = math.floor((b_color_deg - a_color_deg) / scale * delta_b + a_color_deg)

                        if color_deg > 250:
                            img.itemset((y, x, 1), 500 - color_deg)
                            img.itemset((y, x, 0), color_deg - 250)
                        else:
                            img.itemset((y, x, 1), color_deg)
                            img.itemset((y, x, 2), 250 - color_deg)
                    ydraw += scale
                ydraw = 0
            xdraw += scale
        ax.imshow(img)
    @staticmethod
    def tempmap_r(z, ax, t):
        global img
        if z % graphs.scale == 0:
            if z < z_part*10/2:
                z = z + 1
            else:
                z = z - 1
        r_values = img[z].tolist()[::-1]
        whichr = 0.009
        for r_every in r_values:
            rgb_values = list(map(lambda divide255: float(divide255)/255.0, r_every))
            if whichr == 0.009 - 0.009/len(r_values):
                circles = plt.Circle((0, 0), whichr, color=([*rgb_values]), label="T_min={}".format(round(np.min(T[t][math.floor(z/10)]), 4)))
            elif whichr >= 0.009/len(r_values) and whichr <= 0.009*2/len(r_values):
                circles = plt.Circle((0, 0), whichr, color=([*rgb_values]), label="T_max={}".format(round(np.max(T[t][math.floor(z/10)]), 4)))
            else:
                circles = plt.Circle((0, 0), whichr, color=([*rgb_values]))
            ax.add_patch(circles)
            whichr -= 0.009/len(r_values)
        ax.set_title("Temperature Change by Radius in Z={}mm".format(z/10))
        ax.set_xlim([-0.009, 0.009])
        ax.set_ylim([-0.009, 0.009])
        ax.set_xlabel("Radius [m]")
        ax.set_ylabel("Radius [m]")
        ax.legend()
    @staticmethod
    def tempgraph_r(ax, entered_t, entered_z):
        ax.plot(list(map(lambda zx: zx * dr, list(range(0, (r_part + 1))))), T[entered_t][math.floor(entered_z / 10)], label="Temperatures at z = {}mm".format(math.floor(entered_z/10)))
        ax.plot(list(map(lambda zx: zx * dr, list(range(0, (r_part + 1))))), T[entered_t][math.floor(z_part / 2)], label="Center")
        ax.plot(list(map(lambda zx: zx * dr, list(range(0, (r_part + 1))))), T[entered_t][math.floor(z_part)], label="Bottom")
        ax.plot(list(map(lambda zx: zx * dr, list(range(0, (r_part + 1))))), T[entered_t][0], label="Top")
        ax.legend()
        ax.set_title("Temperatures at different z [{}s]".format(entered_t*dt))
        ax.set_xlabel("Radius [m]")
        ax.set_ylabel("Temperature [K]")
    @staticmethod
    def h_plot(ax):
        global h_nearused, h_topused, h_bottomused, old_time
        if old_time == -1:
            h_nearused = [sum(h_mean) / len(h_mean) for h_mean in h_nearused]
            h_topused = [sum(h_mean) / len(h_mean) for h_mean in h_topused]
            h_bottomused = [sum(h_mean) / len(h_mean) for h_mean in h_bottomused]
        ax.set_title("Average Convection Coefficient of Different Surfaces\n Varying with Time")
        x_axis_hgraph = [timesteptoflowtime*max_flowtime/t_part for timesteptoflowtime in list(range(0, len(h_nearused)))]
        ax.plot(x_axis_hgraph, h_topused, label="h_top")
        ax.plot(x_axis_hgraph, h_nearused, label="h_near")
        ax.plot(x_axis_hgraph, h_bottomused, label="h_bottom")
        ax.set_ylabel("Convenction Coefficient (h) [W/m^2K]")
        ax.set_xlabel("Time [s]")
        ax.legend()
    @staticmethod
    def T_minmaxmean__pertime(ax, entered_t):
        # x:Iterate time, y:Temperatures, legends:Min,Max,Mean
        ax.set_title("Temperature Change Until {}s [K]".format(entered_t*dt))
        index_T = [scale_t*dt for scale_t in range(0, entered_t+1)]
        ax.plot(index_T, [np.max(T[time]) for time in range(0, entered_t+1)], label="T_max", color="r")
        ax.plot(index_T, [np.min(T[time]) for time in range(0, entered_t+1)], label="T_min", color="b")

        T_mean = []
        for timemean in range(0, entered_t+1):
            T_everyrmean = 0
            for rmean in range(0, r_part+1):
                if rmean == 0:
                    T_everyrmean += np.mean(T[timemean][::][rmean]) * math.pi * ((dr/2)**2)
                elif rmean == r_part:
                    T_everyrmean += np.mean(T[timemean][::][rmean]) * math.pi * (((rmean)*dr)**2 - ((rmean - 0.5)*dr)**2)
                else:
                    T_everyrmean += np.mean(T[timemean][::][rmean]) * math.pi * (((rmean + 0.5)*dr)**2 - ((rmean - 0.5)*dr)**2)

            T_mean.append(T_everyrmean / (((r_part*dr)**2)*math.pi))

        ax.plot(index_T, [T_mean[time] for time in range(0, entered_t+1)], label="T_mean", color="g")

        #ansys_sure, ansys_max, ansys_min = translater_outfile.Main()
        #ax.plot(ansys_sure, ansys_max, color="m", label="Ansys_max")
        #ax.plot(ansys_sure, ansys_min, color="c", label="Ansys_min")
        #ax.plot(index_T, [-0.00004215963581069570*((time)**2) + 0.04599181744904970000*((time)**1) + 295.49936811233900000000 for time in index_T], label="Experiment", color="black")

        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Temperature [K]")
        ax.legend()


    @staticmethod
    def SpecGraphs(entered_t=5000, entered_z=328, newtime=True):
        global img
        my_dpi = 81.59
        fig = plt.figure(figsize=(1920/my_dpi, 1080/my_dpi), dpi=my_dpi)
        axes = fig.add_gridspec(2, 3)
        ax = [fig.add_subplot(axes[0,0]),
              fig.add_subplot(axes[0,1]),
              fig.add_subplot(axes[0,2]),
              fig.add_subplot(axes[1,0]),
              fig.add_subplot(axes[1,1])]

        fig.canvas.set_window_title("RESULTS OF ANALYSIS")
        graphs.h_plot(ax[0])
        graphs.T_minmaxmean__pertime(ax[1], entered_t)
        graphs.tempgraph_r(ax[2], entered_t, entered_z)
        #zr GRAPH START
        if newtime:
            graphs.tempmap_zr(entered_t, ax[3])
        else:
            ax[3].imshow(img)
        T_minpatch_zr = mpatches.Patch(color="blue", label="T_min={}".format(round(np.min(T[entered_t]), 4)))
        T_maxpatch_zr = mpatches.Patch(color="red", label="T_max={}".format(round(np.max(T[entered_t]), 4)))
        ax[3].legend(handles=[T_minpatch_zr, T_maxpatch_zr], loc="center left", bbox_to_anchor=(1.05, 0.5),
          fancybox=True)
        #zr GRAPH END
        graphs.tempmap_r(entered_z, ax[4], entered_t)

        plt.show()

solvingParams()
BattProp()
dataSavers()
material_air = materials(T_initial)
fou_max = 0
nodeFormulas.convection_type = "_natural" #"_natural", "_forced", "_constant" NOT YET WORKING, ONLY NATURAL

if input("{}s will be divided into {}s intervals, some necessary values for calculation: dz={}m, dr={}m, Fou_z={}. Do you want to keep continuing?[y/n]\n"
                 .format(max_flowtime, dt, dz, dr, alfa*k_find.k_z*dt/(dz**2))).lower() == "n":
    exit()
while t*dt < max_flowtime:
        for j in range(0, z_part+1):
            for i in range(0, r_part+1):
                if (j == 0 or j == z_part) and i == r_part:
                    nodeFormulas.T_convcorner(t, j, i)
                elif i == r_part:
                    nodeFormulas.T_convr(t, j, i)
                elif j == z_part or j == 0:
                    nodeFormulas.T_convz(t, j, i)
                else:
                    nodeFormulas.T_condall(t, j, i)

        t = t + 1
        print("%",round(t/t_part*100, 2))
old_time = -1
while True:
    input_values = input("Entering format of values: [t z], (t:sec, z:meter)").split(" ")
    if len(input_values) == 2:
        try:
            selected_t, selected_z = float(input_values[0]), float(input_values[1])
        except:
            print("Please enter a valid number!")
            continue
    else:
        if input_values[0] == "c":
            break
        else:
            print("You entered non-valid number!")
            continue
    if float(selected_t) <= max_flowtime and float(selected_z) <= 0.065:
        selected_t = int(math.floor(selected_t / dt))  # Seçilen saniyedeki süre süre partına dönüştürülür. Integer dönmesi için yuvarlama yapılır.
        selected_z = int(math.floor(selected_z / dz * graphs.scale))

        graphs.SpecGraphs(selected_t, selected_z, newtime=not (old_time == selected_t))
    else:
        print("Please enter the values: t<{} and z<=0.065".format(max_flowtime))
        continue
    old_time = selected_t
