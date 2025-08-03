using GLMakie
using Colors
using Makie
using LinearAlgebra
using DifferentialEquations
using OrdinaryDiffEq
using NativeFileDialog

#the full function to launch the app
function heat_simulation_app()
    #function that runs the calculations for heat
    function run_calculations(args)
        L, N, α, ρ, c, k, T, steps, solver, profile, base_temp, amp, lbc, lbc_val, rbc, rbc_val = args

        #disecretizes time and length
        dx = L / (N-1)
        dt = T / steps
        x = range(0, L, length=N)
        tspan = (0.0, T)

        #recalculates alpha even though user has provided it just in case
        α = k / (ρ * c)

        #generates the initial conditions based on users input
        u0 = begin
            if profile == :constant
                fill(base_temp, N)
            elseif profile == :sine
                base_temp .+ amp .* sin.(π .* x ./ L)
            elseif profile == :gaussian
                μ, σ = L/2, L/6
                base_temp .+ amp .* exp.(-((x .- μ).^2) ./ (2σ^2))
            elseif profile == :linear
                collect(LinRange(base_temp, base_temp + amp, N))
            end
        end

        #constructs the fininite difference Laplacian matrix
        D = zeros(Float64, N, N)
        for i in 2:N-1
            D[i, i-1] = 1.0
            D[i, i]   = -2.0
            D[i, i+1] = 1.0
        end
        D ./= dx^2  
        
        #sets the boundary conditions based on users input
        bc_type = (lbc == "Dirichlet" ? :dirichlet : :neumann, rbc == "Dirichlet" ? :dirichlet : :neumann)
        bc_vals = (lbc_val, rbc_val)
        if bc_type[1] == :dirichlet
            u0[1] = bc_vals[1]
        end
        if bc_type[2] == :dirichlet
            u0[end] = bc_vals[2]
        end

        #defines the right hand side of the 1d heat equation
        function rhs!(du, u, p, t)
            D, bc_type, bc_val, α = p
            du .= α .* (D * u)

            #defines left boundary condition
            if bc_type[1] == :dirichlet
                du[1] = 0.0
            elseif bc_type[1] == :neumann
                du[1] = α * (u[2] - u[1] - dx * bc_val[1]) / dx^2
            end

            #defines right boundary condition
            if bc_type[2] == :dirichlet
                du[end] = 0.0
            elseif bc_type[2] == :neumann
                du[end] = α * (u[end-1] - u[end] + dx * bc_val[2]) / dx^2
            end
        end

        #creates and solves the heat ode 
        parameters = (D, bc_type, bc_vals, α)
        problem = ODEProblem(rhs!, u0, tspan, parameters)
        solution = solve(problem, solver, dt=dt, saveat=dt, adaptive=false)
        return solution, x
    end

    #creates the screen and the figure
    screen = GLMakie.Screen()
    GLMakie.set_title!(screen, "Heat Rod Simulation - Control")
    fig = Figure(size = (1200, 700), backgroundcolor = colorant"#F0F0F0", font="Times New Roman")
    display(screen, fig)

    #customizes the grid layout of the figure
    for r in 1:10, c in 1:4
        fig[r, c] = Label(fig, "")
    end
    for col in 1:4
        colsize!(fig.layout, col, Relative(0.25))
    end
    rowsize!(fig.layout, 1, Relative(0.15))
    for r in 2:9
        rowsize!(fig.layout, r, Relative(0.0875))
    end
    rowsize!(fig.layout, 10, Relative(0.15))

    #creates a custom theme for the widgets
    custom_theme = Theme(
        Button = (fontsize = 32, buttoncolor = colorant"#F66668", labelcolor = colorant"#8a2929", strokewidth = 0, cornerradius = 100, buttoncolor_active = colorant"#ffffff", buttoncolor_hover =  colorant"#e76163", labelcolor_hover = colorant"#8a2929", labelcolor_active = colorant"#F0F0F0", justification = :center, font="Times New Roman", valign = :bottom, halign = :center),
        Label = (fontsize = 26, color = colorant"#9B3737", justification = :center, font="Times New Roman"),
        Textbox = (fontsize = 26, textcolor = colorant"#9B3737", justification = :left, borderwidth = 0, textcolor_placeholder = colorant"#9B3737", placeholder = "Enter..", cursorcolor = colorant"#9B3737", font="Times New Roman"),
        Menu = (fontsize = 26, textcolor = colorant"#9B3737", font = "Times New Roman", cell_color_active = colorant"#F0F0F0", cell_color_inactive_even = colorant"#F0F0F0", cell_color_inactive_odd = colorant"#F0F0F0", cell_color_hover = colorant"#F66668", strokewidth = 0, direction = :down)
    )
    set_theme!(custom_theme)

    #creates the main label at the top
    main_label_text = Observable("1D Heat Simulation - Enter Parameters")
    main_label = Label(fig, main_label_text, halign = :center, fontsize = 35, valign = :top)
    fig[1, :] = main_label

    #defines the paramaters that the user must input
    param_names = ["Rod Length (m):", "Number of Points (N):", "Thermal Diffusivity (α):", "Density (ρ):", "Specific Heat Capacity (c):", "Thermal Conductivity (k):", "Total Time (s):", "Time Steps:",
        "Method:", "Temp Profile:", "Base Temperature (TK):", "Amplitude (if needed, K):", "Left BC Type:", "Left BC Value:", "Right BC Type:", "Right BC Value:"
    ]

    #defines the solvers that the user can select
    solver_options = Dict("RK4 (Classic)" => RK4(), "Implicit Euler" => ImplicitEuler(), "Trapezoid Rule" => Trapezoid(),
        "KenCarp4 (IMEX)" => KenCarp4(), "Euler (RK1)" => Euler(), "Heun (RK2)" => Heun(), "AB4 (Multistep)" => AB4(),
        "Leapfrog" => PseudoVerletLeapfrog(), "RK9 (Verner)" => Vern9(), "RK8 (DOP)" => DP8()
    )

    #defines the temperature distributions the user can select
    profile_options = Dict("Constant" => :constant, "Sine Wave" => :sine, "Gaussian Pulse" => :gaussian, "Linear Gradient" => :linear)

    #creates the labels that holds the parameter names and the entries to hold them
    entries = []
    global label_idx = 1
    for r in 2:9
        #creates the labels
        fig[r, 1] = Label(fig, param_names[label_idx], halign = :right, justification = :right)
        global label_idx += 1
        fig[r, 3] = Label(fig, param_names[label_idx], halign = :right, justification = :right)
        global label_idx += 1

        #creates the entries
        box1 = Textbox(fig, halign = :left)
        fig[r, 2] = box1
        push!(entries, box1)
        box2 = Textbox(fig, halign = :left)
        fig[r, 4] = box2
        push!(entries, box2)
    end

    #defines the preset materials
    materials = Dict(
        "Aluminum" => (α = 9.7e-5, ρ = 2700.0, c = 897.0, k = 237.0),
        "Copper" => (α = 1.11e-4, ρ = 8960.0, c = 385.0, k = 401.0),
        "Steel" => (α = 1.172e-5,ρ = 7850.0, c = 490.0, k = 43.0)
    )

    #creates the buttons for the presets
    al_btn = Button(fig, label="      Aluminum      ")
    fig[10, 1] = al_btn
    cpr_btn = Button(fig, label="         Copper         ")
    fig[10, 2 ] = cpr_btn
    stl_btn = Button(fig, label="          Steel          ")
    fig[10, 3] = stl_btn

    #function that applies the presets to the correct values
    function apply_material!(entry_vec, material_name)
        mat = materials[material_name]
        Makie.set!(entry_vec[3], string(mat.α))
        Makie.set!(entry_vec[4], string(mat.ρ))
        Makie.set!(entry_vec[5], string(mat.c))
        Makie.set!(entry_vec[6], string(mat.k))
    end

    #connects the buttons to the function
    on(al_btn.clicks) do _
        apply_material!(entries, "Aluminum")
    end

    on(cpr_btn.clicks) do _
        apply_material!(entries, "Copper")
    end

    on(stl_btn.clicks) do _
        apply_material!(entries, "Steel")
    end

    #creates the two dropdown menus for the solver and the temperature profile 
    solver_menu = Menu(fig[6, 2], options = collect(keys(solver_options)), halign = :left, default = first(collect(keys(solver_options))))
    profile_menu = Menu(fig[6, 4], options = collect(keys(profile_options)), halign = :left, default = first(collect(keys(profile_options))))

    #craetes the two dropdown menus for the boundary condition types
    left_bc_menu = Menu(fig[8, 2], options=["Dirichlet", "Nuemann"], halign = :left, default = "Dirichlet")
    right_bc_menu = Menu(fig[9, 2], options=["Dirichlet", "Nuemann"], halign = :left, default = "Dirichlet")

    #the function that confirms that all the inputs are valid
    function begin_calculations!(entries, solver_menu, profile_menu, rbc_menu, lbc_menu)
        args = []
        try
            #makes sure that all the values are entered in the correct 3
            length = parse(Float64, entries[1].displayed_string[])
            points = parse(Int, entries[2].displayed_string[])
            α = parse(Float64, entries[3].displayed_string[])
            ρ = parse(Float64, entries[4].displayed_string[])
            c = parse(Float64, entries[5].displayed_string[])
            k = parse(Float64, entries[6].displayed_string[])
            total_time = parse(Float64, entries[7].displayed_string[])
            time_steps = parse(Int, entries[8].displayed_string[])
            base_temp = parse(Float64, entries[11].displayed_string[])
            amplitude = parse(Float64, entries[12].displayed_string[])
            lbc_val = parse(Float64, entries[14].displayed_string[])
            rbc_val = parse(Float64, entries[16].displayed_string[])

            #checks for nonzero inputs
            if length <= 0 || α <= 0 || ρ <= 0 || c <= 0 || k <= 0 || base_temp < 0 || total_time <= 0 || amplitude < 0
                main_label_text[] = "1D Heat Simulation - Parameters Must Be Positive"
                return nothing
            end

            #checks for error in time step and points
            if points < 3 || time_steps < 3
                main_label_text[] = "1D Heat Simulation - # of Points/Steps Must Be > 2"
                return nothing
            end

            #checks to make sure magnitude is valid
            if amplitude > base_temp
                main_label_text[] = "1D Heat Simulation - Amplitude Cannot Be > Base Temp"
                return nothing
            end

            #after all checking, gets input from the dropdown menus
            solver = solver_options[solver_menu.selection[]]
            profile = profile_options[profile_menu.selection[]]
            rbc = rbc_menu.selection[]
            lbc = lbc_menu.selection[]

            #checks again for errors in boundar conditions
            if (lbc == "Dirichlet" && lbc_val < 0) || (rbc == "Dirichlet" && rbc_val < 0)
                main_label_text[] = "1D Heat Simulation - Dirichlet Conditions Must Be >= 0"
                return nothing
            end

            #defines the final arguments now
            args = [length, points, α, ρ, c, k, total_time, time_steps, solver, profile, base_temp, amplitude, lbc, lbc_val, rbc, rbc_val]
        catch
            #displays an error messgae if invalid
            main_label_text[] = "1D Heat Simulation - Please Ensure Input Are Numbers"
            return nothing
        end

        #calls the functions for calculations
        main_label_text[] = "1D Heat Simulation - Running Simulation"
        solved, x = [], []
        try
            solved, x = run_calculations(args)
        catch e
            main_label_text[] = "1D Heat Simulation - Simulation Failed, Please Try Again"
            print(e)
            return nothing
        end

        #sets up the data
        frames = solved.u
        times = solved.t

        #clears the old fig and then create a new one
        GLMakie.closeall()
        results_screen = GLMakie.Screen()
        GLMakie.set_title!(results_screen, "Heat Rod Simulation - Results")
        results_fig = Figure(size = (1200, 800), backgroundcolor = colorant"#F0F0F0", font="Times New Roman")
        display(results_screen, results_fig)

        #customizes the grid layout of the figure
        for r in 1:12, c in 1:3
            results_fig[r, c] = Label(results_fig, "")
        end
        for col in 1:3
            colsize!(results_fig.layout, col, Relative(0.333333))
        end
        for r in 1:12
            rowsize!(results_fig.layout, r, Relative(0.0825))
        end

        #creates an observable to hold the data of the slider and graph
        temperature_data = Observable(frames[1])

        #creates the the graph that holds the value of temperature vs position at a time
        ax = ax = Axis(results_fig[2:10, :], xlabel="Rod Position (x)", ylabel="Temperature (K)", title="Time = $(round(times[1]; digits=4)) s", xlabelsize=26, ylabelsize=26, titlesize=28, xticklabelsize=22, yticklabelsize=22,
            xgridvisible=false, ygridvisible=false, xminorgridvisible=false, yminorgridvisible=false, bottomspinecolor=:black, topspinecolor=:black, leftspinecolor=:black, rightspinecolor=:black, xtickcolor=:black, 
            ytickcolor=:black, xlabelcolor=colorant"#9B3737", ylabelcolor=colorant"#9B3737", titlecolor=colorant"#9B3737",  xticklabelcolor=colorant"#9B3737", yticklabelcolor=colorant"#9B3737", 
            titlefont = "Times New Roman", xticklabelfont = "Times New Roman", yticklabelfont = "Times New Roman", ylabelfont = "Times New Roman", xlabelfont = "Times New Roman"
        )

        #plots sets up the range of the slider
        lines!(ax, x, temperature_data, linewidth = 9, color = temperature_data, colormap = :reds)
        time_slider_range = LinRange(times[1], times[end], length(times))
        time_slider = Slider(results_fig[11, :], range = time_slider_range, startvalue = times[1], color_active = colorant"#8a2929", color_active_dimmed = colorant"#FF9E9E", color_inactive = colorant"#F66668", linewidth = 10, height = 30)

        #sets the limit of the graph
        alltemps = reduce(vcat, frames)
        ymin = 0.9725 * minimum(alltemps)
        ymax = 1.028277635 * maximum(alltemps)
        ylims!(ax, ymin, ymax)

        #connects the slider value to the results
        on(time_slider.value) do t_val
            idx = searchsortedfirst(times, t_val)

            #renames the title based on time and plots correct values
            if idx <= length(times) && isapprox(times[idx], t_val; atol=1e-8)
                temperature_data[] = frames[idx]
                ax.title = "Time = $(round(times[idx]; digits=4)) s"
            end
        end

        #updates the main label to a new one
        main_label_text[] = "1D Heat Simulation - Adjust Slider to View Results"
        main_results_label = Label(results_fig, main_label_text, halign = :center, fontsize = 33, valign = :top)
        results_fig[1, :] = main_results_label

        #creates the button to save the simulation results
        export_button = Button(results_fig, label = "      Export Results      ")
        results_fig[12, 1] = export_button
        on(export_button.clicks) do _
            @async begin
                path = save_file(filterlist = ".csv;.txt")
                if path !== nothing
                    open(path, "w") do io
                        for (i, u) in enumerate(frames)
                            println(io, "Time=$(times[i]); Temps=", join(u, ","))
                        end
                    end
                end
            end
        end

        #creates the button to exit the program
        exit_button = Button(results_fig, label = "            Exit            ")
        results_fig[12, 2] = exit_button
        on(exit_button.clicks) do _
            exit()
        end

        #creates the button to restart the program
        restart_button = Button(results_fig, label = "          Restart           ")
        results_fig[12, 3] = restart_button
        on(restart_button.clicks) do _
            GLMakie.closeall()
            heat_simulation_app()
        end
    end

    #creates the button that calls the function to begin the calculation
    start_btn = Button(fig, label= "  Begin Calculations  ")
    fig[10, 4] = start_btn
    on(start_btn.clicks) do _
        begin_calculations!(entries, solver_menu, profile_menu, right_bc_menu, left_bc_menu)
    end
end

#calls the function to begin the program
heat_simulation_app()
