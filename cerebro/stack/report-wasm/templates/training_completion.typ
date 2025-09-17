#let recipient = "{{ recipient }}"
#let course = "{{ course }}"
#let date = "{{ date }}"
#let dataset = "{{ dataset }}"
#let logo_path = "logo.png"

#set page(width: 297mm, height: 210mm, margin: 18mm) // A4 landscape
#set heading(numbering: none)
#set par(justify: false)

#let rule = block(height: 1pt, fill: rgb(0,0,0,70%))

#align(center)[
    #image(logo_path, width: 76mm)
    #v(1mm)


    // Title
    #text(size: 42pt, weight: "bold")[Certificate of Completion]
    #v(2mm)
    #text(size: 14pt)[This certifies that]
    #v(2mm)
    #text(size: 28pt, weight: "semibold")[#recipient]
    #v(2mm)
    #text(size: 14pt)[has successfully completed]
    #v(2mm)
    #text(size: 20pt)[#course]


    #text(size: 11pt)[#text(size: 10pt, weight: "bold")[Date:] #date]
    #v(0.5mm)
    #text(size: 11pt)[#text(size: 10pt, weight: "bold")[Dataset:] #dataset]

]