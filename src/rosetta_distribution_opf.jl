module rosetta_distribution_opf

    import JuMP
    import InfrastructureModels
    import InfrastructureModels: optimize_model!, @im_fields, nw_id_default, ismultinetwork, update_data!
    import PowerModelsDistribution
    import LinearAlgebra: diag, diagm
    import LaTeXStrings: latexstring, text

    const IM = InfrastructureModels
    const PMD = PowerModelsDistribution

    using GraphPlot
    using Graphs
    const _GP = GraphPlot
    const _G = Graphs

    import Plots
    export build_mc_opf

    include("./3wire/ACP.jl")
    include("./3wire/ACR.jl")

    include("./4wire_IVR/PMD/base.jl")
    include("./4wire_IVR/PMD/variables.jl")
    include("./4wire_IVR/PMD/constraint_template.jl")
    include("./4wire_IVR/PMD/constraint.jl")
    include("./4wire_IVR/PMD/ivr.jl")

    include("./4wire_IVR/RPMD/IVR_EN.jl")
    include("./4wire_IVR/RPMD/IVR_EN_vectorized.jl")

    include("./util/graph.jl")
    include("./util/helper_functions.jl")
    include("./util/solution.jl")
    include("./util/visualisation.jl")
    

    include("./inverters/variable.jl")
    include("./inverters/constraint_template_en.jl")
    include("./inverters/en_ivr.jl")
    include("./inverters/objective.jl")

    include("./inverters/opf.jl")
    include("./inverters/opf_sop.jl")
    include("./inverters/loss_branch_inverter.jl")
    include("./inverters/loss_branch_sop.jl")

    
end # module rosetta_distribution_opf
