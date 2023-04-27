# code for read intances .txt

using JuMP, CPLEX

struct InstanceData
  N::Int
  v
  sigma
  M::Int
  kr
  Q
  s
  rho
  alpha
end

# define instance name
instance="port6.txt"
path="../instances/$(instance)"

#define results file name
#solution="solutions.txt"

# function for read and print instance
function read_data(path)

  file = open(path)
  fileText = read(file, String)
  tokens = split(fileText) 

  aux = 1
  #read the problem dimension
  N = parse(Int,tokens[aux])
  println(N)

  #aux += 1
  #M = parse(Int,tokens[aux])
  M = Int((N*(N-1)/2 + N))
  println(M)

  # rho = 0.0018
  aux += 1
  rho = parse(Float64,tokens[aux])
  println(rho)

  #alpha = 15
  aux += 1
  alpha = parse(Float64,tokens[aux])
  println(alpha)

  v = zeros(Float64,N)  
  for i in 1:N
    aux += 1
    v[i] = parse(Float64,tokens[aux])
  end

  for i in 1:N
    v[i] = 4*v[i]
  end

  #println(v)
    
  sigma = zeros(Float64,N)  
  for i in 1:N
    aux += 1
    sigma[i] = parse(Float64,tokens[aux])
  end
  ##println(sigma)

  kr = zeros(Float64,M)  
  for i in 1:M
    aux += 1
    kr[i] = parse(Float64,tokens[aux])
  end
  ##println(kr)
  
  
  s = zeros(Int64,N)
  s[1] = 0
  for i in 2:N
    s[i] = N + 0.5*((2*N-1)-i+2)*(i-2) # 61?
  end

  Q = zeros(Float64,N,N)  
  for i in 1:N
    for j in i:N
        Q[i,j] = kr[s[i]+j-i+1]*sigma[i]*sigma[j]
        Q[i,j] = 4*Q[i,j]
        Q[j,i] = Q[i,j]
    end 
  end
  
  close(file) # close file instances
  
  inst = InstanceData(N,v,sigma,M,kr,Q,s,rho,alpha)
  
  return inst
  
end

# function to model the problem 
function main(inst::InstanceData)
  maxtime = 3600
  tolgap = 0.000001

  model = Model(CPLEX.Optimizer)
  
  set_optimizer_attribute(model,"CPX_PARAM_TILIM",maxtime) # time limit
  set_optimizer_attribute(model,"CPX_PARAM_EPGAP",tolgap) # relative MIP optimality gap
  #set_optimizer_attribute(model,"CPX_PARAM_LPMETHOD ",0) # method used in root node
  #set_optimizer_attribute(model,"CPX_PARAM_NODELIM",maxnodes) # MIP node limit
  #set_optimizer_attribute(model,"CPX_PARAM_THREADS",1) # number of threads 
  set_optimizer_attribute(model, "CPXPARAM_OptimalityTarget", 3)   

  N = inst.N
    
  @variable(model, y[i=1:N], Bin)
  @variable(model, 0<=x[i=1:N]<=1)
  
  @objective(model, Min, sum(x[i]*inst.Q[i,j]*x[j] for i in 1:N, j in 1:N))
  
  @constraint(model, constraint1[t in 1:1], sum(inst.v[i]*x[i] for i in 1:N) >= inst.rho)
  @constraint(model, constraint2[t in 1:1], sum(x[i] for i in 1:N) == 1)
  @constraint(model, constraint3[t in 1:1], sum(y[i] for i in 1:N) >= inst.N - inst.alpha)
  @constraint(model, constraint4[i in 1:N], x[i] + y[i] <= 1)
  
  method = "mip"
  if method == "lp"
    undo_relax = relax_integrality(model)
  end

  #write_to_file(model,"problem.lp")

  JuMP.optimize!(model)
  
  ### solving the problem ###
  status = optimize!(model)
  
  opt = 0
  if termination_status(model) == MOI.OPTIMAL    
    println("status = ", termination_status(model))
    opt = 1
  else
    println("status = ", termination_status(model))
  end
  
    ### get solutions ###
  
  if method == "mip"
    bestbound = objective_bound(model)
    numnodes = node_count(model)
    gap = MOI.get(model, MOI.RelativeGap())
  end
  bestsol = objective_value(model)
  time = solve_time(model)

  ### print results ###
  open("../results/result_$(instance)","a") do f
    if method == "mip"
      write(f,"$(instance);$(N);$(method);$(bestbound);$(bestsol);$(gap);$(time);$(numnodes);$(opt)\n")
    else
      write(f,"$(instance);$(N);$(form);$(method);$(bestsol);$(time)\n")
    end
  end

#  open("sol_$(instance)", "w") do f
#    println(f,"x:")
#    for i in 1:N
#      if x[i] != 0.0
#        write(f,"x[$i] = ", JuMP.value(x[i]))
#      end
#    end

#    println(f,"y:")
#    for i in 1:N
#      if y[i] == 0
#        println(f,"y[$i] = ", JuMP.value(y[i]))
#      end
#    end    
#  end

end

# call function
inst = read_data(path)

main(inst)
