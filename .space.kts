/**
* JetBrains Space Automation
* This Kotlin-script file lets you automate build activities
* For more info, see https://www.jetbrains.com/help/space/automation.html
*/

job("check shinyGatom package") {
    git("shinyCyJS")

    container(displayName = "run tests", image = "ctlab.registry.jetbrains.space/p/shiny-apps/containers/shiny-apps-base") {
    	shellScript {
        	content = """
                R -e 'devtools::install("'$mountDir/work/shinyCyJS/'", upgrade=FALSE, quick=TRUE)'
                R -e "devtools::install_dev_deps(upgrade=FALSE)"
                R -e "devtools::test(stop_on_failure=TRUE)"
            """
        }
    }
    
}

