import random
import numpy as np
import math as m
matrix=np.loadtxt("DAG.txt", dtype='i', delimiter=',')
print(matrix)
(r,c)=matrix.shape  #rows and cols
Tmatrix=np.array(matrix).transpose() #transpose matrix
#print()
#print(Tmatrix)
# for entry tasks finding from DAG---------------------------------------
en_task=[]
for i in range(r):
    count=0
    for j in range(c):
        if Tmatrix[i][j]==0:
            count+=1
    if count==r:
        en_task.append(i)
print("Entry tasks are:",en_task)
#for exit tasks finding from DAG-------------------------------------
ex_task=[]
for i in range(r):
    count=0
    for j in range(c):
        if matrix[i][j]==0:
            count+=1
    if count==r:
        ex_task.append(i)
print("exit_tasks are:",ex_task)
#for predecessor tasks
pred=[]
for i in range(r):
    emp=[]
    for j in range(c):
        if Tmatrix[i][j]>0:
            emp.append(j)
    pred.append(emp)
print("predecessors list for tasks labelled from 0 to "+str(len(pred))+" are:",pred)
for i in range(len(pred)):
    print("predecessor tasks for Task_"+str(i)+" are:",pred[i])
#for finding Successor of tasks
succ=[]
for i in range(r):
    emp=[]
    for j in range(c):
        if matrix[i][j]>0:
            emp.append(j)
    succ.append(emp)
print("successors list for tasks labelled from 0 to "+str(len(succ))+" are:",succ)
for i in range(len(pred)):
    print("successor tasks for Task_"+str(i)+" are:",succ[i])
#length of tasks Given In the DAG
Task_len=[]
DAG_len= np.loadtxt("DAG_len.txt", dtype='i', delimiter=',')
for i in range(r):
    Task_len.append(DAG_len[i])
print(Task_len)
#Generating Tasks in different Levels
su=[]
Lev_tasks=[]
for i in range(len(succ)):
    if i==0:
        su.append([i])
    else:
        emp=[]
        for j in range(len(su[i-1])):
            x=su[i-1][j]
            for k in range(len(succ[x])):
                emp.append(succ[x][k])
        su.append(list(set(emp)))
for k in range(len(su)):
    if len(su[k]) !=0:
        Lev_tasks.append(sorted(su[k]))
print("No of Levels In Given DAG :", len(Lev_tasks))
for k in range(len(Lev_tasks)):
    print("Tasks in Level_",str(k),"are :",Lev_tasks[k])
#for Generating the Population
pop=[]
for i in range(len(Lev_tasks)):
    for k in range(len(Lev_tasks[i])):
        pop.append(Lev_tasks[i][k])
print(pop)
#for calculating pe
def pe(v,ex):
        beta = 0.66  # constant
        IV = 230  # in volts given
        D=0.24
        DV = [0, 0.25, 0.5, 1.0, 1.25]  # Decrement in voltage given
        AFT = []
        EST = []
        ET = []
        VM_assign = []
        avail = [0] * len(sp_dv)
        for j in range(len(v)):
            et_emp = []
            for k in range(len(sp_dv)):
                exec_time = (Task_len[v[j]]) / (sp_dv[k] * (1 - D))  # ET calculation i calculated as per formula given in objectives
                et_emp.append(exec_time)
            if j == 0:
                EST.append(0)
                AFT.append(min(et_emp) + EST[j])
                ET.append(min(et_emp))
                ind = et_emp.index(min(et_emp))
                VM_assign.append(ind)
                avail[ind] = AFT[ind]
            else:
                emp = []
                ele = v[j]
                ind1 = v.index(v[j])
                for k in range(len(pred[j])):
                    emp.append(AFT[pred[j][k]])
                EST.append(max(emp))
                c = 0
                for k in range(len(avail)):
                    if avail[k] <= EST[j] and c == 0:
                        VM_assign.append(k)
                        ET.append(et_emp[k])
                        AFT.append(et_emp[k] + EST[j])
                        avail[k] = AFT[j] - ET[j]
                        c += 1
        mks = AFT[ex_task[0]]  # minized makespan of molecule
        N_mks=mks / sum(AFT)  # normalized makespan
        Econ = 2 * mks * beta * IV * sum(DV)  # energy conserved of each molecule
        N_Econ=Econ/ex
        gamma = 0.7  # any value between 0 and 1 and delta =1-gamma
        delta = 1 - gamma
        PE=((gamma * N_mks) + (delta / (1 + N_Econ)))
        return  PE
#for calculating the Econ ,mks of reuired molecules
def ms(v, ex):
    beta = 0.66  # constant
    IV = 230  # in volts given
    D = 0.24
    DV = [0, 0.25, 0.5, 1.0, 1.25]  # Decrement in voltage given
    AFT = []
    EST = []
    ET = []
    VM_assign = []
    CPmin=[]
    avail = [0] * len(sp_dv)
    for j in range(len(v)):
        et_emp = []
        for k in range(len(sp_dv)):
            exec_time = (Task_len[v[j]]) / (sp_dv[k] * (1 - D))  # ET calculation i calculated as per formula given in objectives
            et_emp.append(exec_time)
        CPmin.append(min(et_emp))
        if j == 0:
            EST.append(0)
            AFT.append(min(et_emp) + EST[j])
            ET.append(min(et_emp))
            ind = et_emp.index(min(et_emp))
            VM_assign.append(ind)
            avail[ind] = AFT[ind]
        else:
            emp = []
            ele = v[j]
            ind1 = v.index(v[j])
            for k in range(len(pred[j])):
                emp.append(AFT[pred[j][k]])
            EST.append(max(emp))
            c = 0
            for k in range(len(avail)):
                if avail[k] <= EST[j] and c == 0:
                    VM_assign.append(k)
                    ET.append(et_emp[k])
                    AFT.append(et_emp[k] + EST[j])
                    avail[k] = AFT[j] - ET[j]
                    c += 1
    mks = AFT[ex_task[0]]  # minized makespan of molecule
    N_mks = mks / sum(AFT)  # normalized makespan
    Econ = 2 * mks * beta * IV * sum(DV)  # energy conserved of each molecule
    N_Econ = Econ / ex
    gamma = 0.7  # any value between 0 and 1 and delta =1-gamma
    delta = 1 - gamma
    PE = ((gamma * N_mks) + (delta / (1 + N_Econ)))
    return mks,Econ,AFT,ET,EST,VM_assign,CPmin
#Algorithm-2 Implemenation -------------------------------
def ON_WALL_PSEUDO_EFFECTIVE(V_i1):
    V_new=V_i1.copy()
    T_i = random.choice(V_new[1:r - 1])
    ind=V_new.index(T_i)
    T_j = pred[T_i][0]          # last predecessor
    T_k = succ[T_i][0]            # first successor
    PT_i=V_new.index(T_j)             #position of T_i and T_j
    PT_j=V_new.index(T_k)
    rpos=random.randint(PT_i+1,PT_j)  #rpos!=PT_i #print(PT_i,PT_j) #print(rpos)---to check use these prints
    T_dash=V_new[rpos]                    #atom at rpos in V_new
    ind1=V_new.index(T_dash)
    V_new[ind]=T_dash            #exchange T_dash with T_i in V_new
    V_new[ind1]=T_i              #exchange T_i with T_dash
    #print(V_new)
    return V_new
#Algorithm-3 Implemenatation --------------------------------
def REPLACE(v,pop_size,V_neww,PE,ex,KE):
    print("Before Modification ")
    print(v)
    V_Best=v[0].copy()
    buf=[] #which stores best molecules
    buf.append(V_Best)
    V_Worst=v[0].copy()
    PE_Best=PE[0]
    PE_Worst=PE[0]
    PE_Vnew=pe(V_neww,ex)
    for i in range(pop_size):
        if PE[i]<=PE_Best:        #print(v[i])#print(PE[i]) use to check
            V_Best=v[i].copy()
            PE_Best=PE[i]
            ind=i
        elif PE[i]>=PE_Worst:
            V_Worst=v[i].copy()   ##print(V_Best,V_Worst) use to check
            PE_Worst=PE[i]
            ind1=i
    if PE_Vnew<=PE_Best:
        V_Best=V_neww.copy()
        v[ind]=V_Best  #replacing best in population
        buf.append(V_Best)
    elif PE_Vnew<=PE_Best+KE:
        V_Best = V_neww.copy()
        v[ind] = V_Best  # replacing best in population
    elif PE_Vnew<=PE_Worst+KE:
        V_Worst=V_neww.copy()
        v[ind1]=V_Worst
    return v
#Algorith-1 Implemenation --------MAIN ALGORITHM------------
def EnergyEfficientWorkFlow(pop,type,intertype,KE_initial,D,sp_dv):
    max_i = 3  # maximum no of iterations we need to define
    it=0
    while it<=max_i:
                                                                       #Generating random populatin
        print("Population after Iteration ",str(it)," ---------------------------------------------------------------------------------------------")
        print()
        c_pop = pop.copy()
        c_pop.remove(ex_task[0])
        beta = 0.66  # constant
        buff=0.3 #taken
        IV = 230  # in volts given
        KE=0.1 #ke initial
        DV = [0, 0.25, 0.5, 1.0, 1.25]  # Decrement in voltage given
        pop_size = 4  # population size
        v = []
        for i in range(pop_size):
            atom = random.choice(c_pop[1:r - 1])
            # print(atom)
            c_pop.remove(atom)
            l_pred = pred[atom][0]  # last predecessor
            f_suc = succ[atom][0]  # first successor
            l = pop.copy()
            z = l[atom:f_suc]
            z.reverse()
            l = l[en_task[0]:atom] + z + l[f_suc:ex_task[0] + 1]
            v.append(l)
        for i in range(len(v)):
            print("Molecule ", str(i), ":", v[i])
        ex = 0
        N_mks = []
        Econ_temp = []
        N_Econ = []
        PE = []
        for i in range(len(v)):
            mol = v[i]
            AFT = []
            EST = []
            ET = []
            VM_assign = []
            avail = [0] * len(sp_dv)
            for j in range(len(v[i])):
                et_emp = []
                for k in range(len(sp_dv)):
                    exec_time = (Task_len[v[i][j]]) / (
                                sp_dv[k] * (1 - D))  # ET calculation i calculated as per formula given in objectives
                    et_emp.append(exec_time)
                if j == 0:
                    EST.append(0)
                    AFT.append(min(et_emp) + EST[j])
                    ET.append(min(et_emp))
                    ind = et_emp.index(min(et_emp))
                    VM_assign.append(ind)
                    avail[ind] = AFT[ind]
                else:
                    emp = []
                    ele = v[i][j]
                    ind1 = v[i].index(v[i][j])
                    for k in range(len(pred[ele])):
                        emp.append(AFT[pred[ele][k]])
                    EST.append(max(emp))
                    c = 0
                    for k in range(len(avail)):
                        if avail[k] <= EST[j] and c == 0:
                            VM_assign.append(k)
                            ET.append(et_emp[k])
                            AFT.append(et_emp[k] + EST[j])
                            avail[k] = AFT[j] - ET[j]
                            c += 1
            mks = AFT[ex_task[0]]  # minized makespan of molecule
            N_mks.append(mks / sum(AFT))  # normalized makespan
            Econ = 2 * mks * beta * IV * sum(DV)  # energy conserved of each molecule
            Econ_temp.append(Econ)
            ex += Econ
            print()
            print("For Molecule ", str(i), "--------")
            print()
            print("Vm assigned:", VM_assign)
            print("AFT :", AFT)
            print("EST :", EST)
            print("ET :", ET)
            print("Minimum Makespan of Each molecule :", mks)
            print("Maximum Energy Conserved Of Each Molecule :", Econ)
        for k in range(len(Econ_temp)):
            N_Econ.append(Econ_temp[k] / ex)  # normalized Energy conserved values
        gamma = 0.7  # any value between 0 and 1 and delta =1-gamma
        delta = 1 - gamma
        for k in range(len(Econ_temp)):
            PE.append((gamma * N_mks[k]) + (delta / (1 + N_Econ[k])))
            print("Potential Energy of Molecule ", str(k), " :", PE[k])
        V_i1=[]
        V_i2=[]
        best=[]
        best_pe=[]
        for i in v:
            ind=v.index(i)
            r1 =random.uniform(0.1, 0.8) #choosing random number r1
            if r1<=type:
                v_OWIC = i.copy()# molecule to perform on wall ineffective collision
                atom1 = random.choice(v_OWIC[1:r - 1])
                #print(atom1)
                la_pred = pred[atom1][0]
                z1 = v_OWIC[la_pred + 1:atom1 + 1]
                #print(z1)
                z1.reverse()
                v_OWIC = v_OWIC[en_task[0]:la_pred+1] + z1 + v_OWIC[atom1 + 1:ex_task[0] + 1]  #molecules after onwall ineffective collision
                #print(v_OWIC)
                pe_of_vOWIC=pe(v_OWIC,ex) #potential enrgy of above molecule
                # print(pe_of_vOWIC)
                if pe_of_vOWIC<=PE[ind]+KE: #condition check for On wall ineffective collision
                    v_i_1=v_OWIC.copy()
                pe_vi_1=pe(v_i_1,ex)
                #print(v_i_1)
                # decomposition on v_i1
                atom2=random.choice(v_i_1[1:r - 1])
                la1_pred = pred[atom2][0]
                z2= v_i_1[la1_pred + 1:atom2 + 1]
                z2.reverse()
                v_i_2 = v_i_1[en_task[0]:la1_pred + 1] + z2 + v_i_1[atom2 + 1:ex_task[0] + 1] #after decomposition
                pe_vi_2 = pe(v_i_2, ex)
                # print(v_i_2)
                #print(pe_vi_2)
                if pe_vi_1 +pe_vi_2<=PE[ind]+KE+buff:
                    V_i1=v_i_2.copy()
                    print("Molecule after performing on wall ineffective and decomposition V_i","1 :",V_i1)
            else:
                v_k=[] #selecting another molecule from parent
                v_i=[]
                v_i=i.copy()
                s= random.choice(v)
                v_k=s.copy()
                pe_of_vi=pe(v_i,ex)
                pe_of_vk=pe(v_k,ex)
                #print(v_i,v_k)
                r2 = random.uniform(0.1, 0.8) #choosing random number r2
                if r2<=intertype:
                    #performing intermolecular ineffective collsion
                    cr_point=random.choice(v_k[1:r - 1]) #randomly choosing crossover point
                    new_list=v_i[cr_point:]
                    new_list_1=v_k[cr_point:]
                    for k in range(len(new_list)):
                        if new_list[k] not in new_list_1:
                            t=new_list[k]
                            new_list[k]=new_list_1[k]
                            new_list_1[k]=t
                    v_i_dash=v_i[en_task[0]:cr_point]+new_list_1  # 1- point crossover operation
                    v_k_dash=v_k[en_task[0]:cr_point]+new_list    #print(v_i_dash) #print(v_k_dash)#print("aftercross:") #use print to check
                    pe_of_vi_dash=pe(v_i_dash,ex)
                    pe_of_vk_dash=pe(v_k_dash,ex)   #print(pe_of_vi,pe_of_vk) use it to check
                    if pe_of_vi_dash <=pe_of_vk_dash:
                        if pe_of_vk_dash+pe_of_vi_dash <=pe_of_vk+pe_of_vi+KE+KE:    #condtion check
                            V_i1=v_i_dash.copy()
                    elif pe_of_vk_dash <=pe_of_vi_dash:
                        if pe_of_vk_dash+pe_of_vi_dash <=pe_of_vk+pe_of_vi+KE+KE:
                            V_i1=v_k_dash.copy()
                    print("Molecule after intermolecular Ineffective collision :",V_i1)
                else:
                    #systhesis operation here we collide two molecules using crossover
                    cr_point = random.choice(v_k[1:r - 1])
                    new_list=v_i[cr_point:]
                    new_list_1=v_k[cr_point:]
                    for k in range(len(new_list)):
                        if new_list[k] not in new_list_1:
                            t = new_list[k]
                            new_list[k] = new_list_1[k]
                            new_list_1[k] = t
                    v_i_syn = v_i[en_task[0]:cr_point] + new_list_1
                    v_k_syn = v_k[en_task[0]:cr_point] + new_list
                    pe_of_vi_syn = pe(v_i_syn, ex)
                    pe_of_vk_syn = pe(v_k_syn, ex)
                    #nowcondition check we select lower pe molecule
                    if pe_of_vi_syn <=pe_of_vk_syn :
                        if pe_of_vi_syn<=pe_of_vk+pe_of_vi+pe_of_vk_syn+KE+KE: #conditoin check
                            V_i1=v_i_syn.copy()
                    elif pe_of_vk_syn<pe_of_vi_syn :
                        if pe_of_vk_syn <= pe_of_vk + pe_of_vi + pe_of_vi_syn + KE + KE:  # conditoin check
                              V_i1 = v_k_syn.copy()
                    print("Molecule After Sysnthesis operation :",V_i1)
            V_i2=ON_WALL_PSEUDO_EFFECTIVE(V_i1) #calling Algo-2
            PE_V_i1=pe(V_i1,ex)             #PE of V_i1 and V_i2
            PE_V_i2=pe(V_i2,ex)
            if PE_V_i1<=PE_V_i2:
                V_i=V_i1
            else:
                V_i=V_i2
            Modified_pop=REPLACE(v,pop_size,V_i,PE,ex,KE)  #here v is my population and calling Algorithm-3
            print("population After Modified by performing all operations")
            print(Modified_pop)
            Final_pe=[]
            for k in range(len(Modified_pop)):
                Final_pe.append(pe(Modified_pop[k],ex))
                print("PE of Molecule_",str(k),"in Modified ",Final_pe[k])
            ind=Final_pe.index(min(Final_pe))
            Best_molecule=Modified_pop[ind]
            best.append(Best_molecule)
            best_pe.append(Final_pe[ind])
            print("Molecule_", str(ind), " is the best")
            print("Molecule is :", Best_molecule)
        min_ind=best_pe.index(min(best_pe))
        print("Final the best molecule after all iterations : ",best[min_ind])
        print("Its fitness value is :",pe(best[min_ind],ex))
        mks, Econ, AFT, ET, EST, VM_assign,CPmin=ms(best[min_ind],ex)
        # Factors calculation part for best molecule
        #1.makespan
        print("makespan of Best molecule :",mks)
        #2.Energy conserved
        print("Energy conserved of Best Molecule :",Econ)
        #3.speedup
        EET = sum(ET)
        speedup = EET / mks
        print("Speedup is :",round(speedup))
        #4.Vm Utilization
        Vm_util = []
        AWT = []  # actual Working Time
        for i in range(len(sp_dv)):
            ss = []
            for j in range(r):
                if i == VM_assign[j]:
                    ss.append(ET[j])
                else:
                    ss.append(0)
            #print("Execution time of Tasks Executing in Vm " + str(i) + " :", ss)
            if ss == []:
                ss.append(0)
                AWT.append(sum(ss)) #exec time sum of vm
            else:
                Vm_util.append(max(ss) / mks)
                AWT.append(sum(ss))
        print("Actual Working Time :", AWT)
        print("VM utilization Individually respectively :", Vm_util)
        Avg_Vm_util = sum(Vm_util) / len(sp_dv)  # average vm utilization
        print("Average cloud Utilization is :",Avg_Vm_util)
        #5.Load Balancing Factor
        import statistics as st
        LB_sum_1 = 0
        for i in range(len(Vm_util)):
            LB_sum_1 += ((Vm_util[i] - Avg_Vm_util) ** 2)
            SD = m.sqrt(LB_sum_1 / len(sp_dv))
            LB_2 = (abs(Avg_Vm_util - SD) / Avg_Vm_util) * 100
        print("LoadBalacing Factor  : ", LB_2, "%")
        #6.SLR
        SLR=mks/sum(CPmin)
        print("SLR :",SLR)
        it+=1
sp_dv=[4,2.5,2]
EnergyEfficientWorkFlow(pop, 0.4, 0.5, 0.1,0.24,sp_dv)









