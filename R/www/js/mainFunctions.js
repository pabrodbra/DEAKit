/*
window.pressed = function(){
    var a = document.getElementById('fileIn');
    if(a.value === "")    {
        noFile.innerHTML = "No folder choosen";
    }
    else    {
        noFile.innerHTML = "";
    }
};

document.getElementById("folderpicker").addEventListener("change", function(event) {
    let output = document.getElementById("listing");
    let files = event.target.files;
    var arr = [];

    for (let i=0; i<files.length; i++) {
        console.log(files[i])
        arr.push(files[i])
    };

    Shiny.onInputChange("mydata", arr);
}, false);


                       # 6 - FIELD OF VIEW THRESHOLD - For FOV plot
                       numericInput("key.of.interest", "Key",
                                    placeholder = "Number of key in name to compare"),
                       helpText("Path of the Location of the RLF File"),
                       
                       # 7 - BINDING DENSITY THRESHOLDS (min and max) - For BD plot
                       numericInput("key.of.interest", "Key",
                                    placeholder = "Number of key in name to compare"),
                       helpText("Path of the Location of the RLF File"),
                       
                       # 8 - P VALUE THRESHOLD - For DEA & Pathway Enrichment
                       numericInput("key.of.interest", "Key",
                                    placeholder = "Number of key in name to compare"),
                       helpText("Path of the Location of the RLF File"),
                       
                       # 9 - LOG FC THRESHOLD - For DEA & Pathway Enrichment
                       numericInput("key.of.interest", "Key",
                                    placeholder = "Number of key in name to compare"),
                       helpText("Path of the Location of the RLF File"),
                       
                       # 10 - Q VALUE THRESHOLD - For Pathway Enrichment
                       numericInput("key.of.interest", "Key",
                                    placeholder = "Number of key in name to compare"),
                       helpText("Path of the Location of the RLF File"),
                       
                       # 11 - SEED - Optional
                       numericInput("key.of.interest", "Key",
                                    placeholder = "Number of key in name to compare"),
                       helpText("Path of the Location of the RLF File"),
                       
                       # 12 - OUTPUT - Optional
                       textInput("rlf.path", "RLF File",
                                 placeholder = "/home/userXY/data/RLF_file.rlf"),
                       helpText("Path of the Location of the RLF File"),

*/