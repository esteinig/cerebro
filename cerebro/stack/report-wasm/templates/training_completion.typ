#let recipient = "{{ recipient }}"
#let course = "{{ course }}"
#let date = "{{ date }}"
#let dataset = "{{ dataset }}"
#let sensitivity = "{{ sensitivity }}"
#let specificity = "{{ specificity }}"
#let logo_path = "logo.png"

#set page(width: 297mm, height: 210mm, margin: 18mm) // A4 landscape
#set heading(numbering: none)
#set par(justify: false)

#align(center)[
    #image(logo_path, width: {{ logo_width }})
    #v(1mm)

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
    #v(1mm)
    #text(size: 11pt)[#text(size: 10pt, weight: "semibold")[Dataset:] #dataset]
    #h(1.5mm)
    #text(size: 11pt)[#text(size: 10pt, weight: "semibold")[Sensitivity:] #sensitivity]
    #h(1.5mm)
    #text(size: 11pt)[#text(size: 10pt, weight: "semibold")[Specificity:] #specificity]

]