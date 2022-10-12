module test_misc

using Test, XtalsPyTools

@testset "pydep functions" begin
    @test isnothing(XtalsPyTools.load_pydep("bogus_python_package"))
    @test ismissing(
        XtalsPyTools.check_pydep(:bogus_python_package => "bogus_python_package")
    )
end

end
