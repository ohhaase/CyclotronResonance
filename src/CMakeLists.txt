# This file is being included in ./CMakeLists.txt

# Add source files to global project
target_sources(${PROJECT_NAME}
    PRIVATE
        # ./src/
        main.cpp
        global_vars.cpp

        # ./src/core/
        core/helper_functions.cpp
        core/parallel_functions.cpp
        core/sim_types.cpp

        # ./src/helper_objects/
        helper_objects/ElectronDistb.cpp
        helper_objects/Histogram.cpp
        helper_objects/Histogram2D.cpp
        helper_objects/NumericalDistb.cpp
        helper_objects/Timer.cpp

        # ./src/science_functions/
        science_functions/beta_calculations.cpp
        science_functions/cross_sections.cpp
        science_functions/sampling_functions.cpp
        science_functions/scattering_functions.cpp

    PUBLIC
        FILE_SET HEADERS
        FILES
            # ./src/
            global_vars.hpp

            # ./src/core/
            core/helper_functions.hpp
            core/parallel_functions.hpp
            core/sim_types.hpp

            # ./src/helper_objects/
            helper_objects/ElectronDistb.hpp
            helper_objects/Histogram.hpp
            helper_objects/Histogram2D.hpp
            helper_objects/NumericalDistb.hpp
            helper_objects/Timer.hpp

            # ./src/science_functions/
            science_functions/beta_calculations.hpp
            science_functions/cross_sections.hpp
            science_functions/sampling_functions.hpp
            science_functions/scattering_functions.hpp
)