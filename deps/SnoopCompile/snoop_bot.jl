using CompileBot

bot_config = BotConfig(
    "Marta";
    os = ["linux", "macos", "windows"],
    else_os = "linux",
    version = [v"1.5.0"],
    else_version = v"1.5.0",
)

snoop_bot(bot_config, "$(@__DIR__)/example_script.jl")
