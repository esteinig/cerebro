
import { FileTag, Role } from "./types";
import type { Cerebro } from "./types";

export function getCssVariableAsHex(variableName: string | null | undefined, theme: string): string | null {
    if (!variableName) return null
    const element = document.querySelector(`[data-theme="${theme}"]`); // Target the themed element
    if (!element) return null;
    const rgbValue = getComputedStyle(element).getPropertyValue(variableName).trim();
    if (!rgbValue) return null;
    const [r, g, b] = rgbValue.split(' ').map(Number);
    return `#${((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1)}`;
}

export const getDateTime = (isoString: string): [string, string] => {
    let data = isoString.split("T");
    return [data[0], data[1].split("Z")[0].substring(0,8)]
}

export const getDateTimeString = (
    isoString: string, 
    time: boolean = true, 
    sep: string = " ", 
    localeString: boolean = false
): string => {
    try {
        const date = new Date(isoString);

        if (localeString) {
            const options: Intl.DateTimeFormatOptions = {
                year: 'numeric',
                month: 'numeric',
                day: 'numeric',
                ...(time && { hour: '2-digit', minute: '2-digit', second: '2-digit' }),
            };
            return date.toLocaleString(undefined, options).replace(',', sep);
        } else {
            const datePart = date.toISOString().split("T")[0];
            const timePart = date.toISOString().split("T")[1].substring(0, 8);
            return time ? `${datePart}${sep}${timePart}` : datePart;
        }
    } catch (error) {
        return isoString;
    }
};

export const getDateTimeStringUtc = (
    utcString: string, 
    time: boolean = true, 
    localeString: boolean = false,
    sep: string = " ", 
): string => {
    const date = new Date(utcString);

    if (localeString) {
        const options: Intl.DateTimeFormatOptions = {
            year: 'numeric',
            month: 'long',
            day: 'numeric',
            ...(time && { 
                hour: '2-digit', 
                minute: '2-digit', 
                second: '2-digit', 
                hour12: false // Use 24-hour format
            }),
        };
        let formattedString = date.toLocaleString(undefined, options);
        
        if (time) {
            formattedString = formattedString.replace(' at ', " @ ");
        }
        return formattedString;
        
    } else {
        const datePart = date.toISOString().split("T")[0];
        const timePart = date.toISOString().split("T")[1].substring(0, 8);
        return time ? `${datePart}${sep}${timePart} UTC` : datePart;
    }
}

export const isWithinTimeLimit = (datetimeStr: string | undefined, minutes: number): boolean => {

    if (datetimeStr === undefined) {
        return false
    }

    // Parse the datetime string to a Date object
    const parsedDate = new Date(datetimeStr);

    // Get the current UTC time
    const currentTime = new Date();

    // Calculate the difference in milliseconds
    const timeDifference = currentTime.getTime() - parsedDate.getTime();

    // Convert the time difference from milliseconds to minutes
    const differenceInMinutes = timeDifference / (1000 * 60);

    // Check if the difference is within the given limit
    return differenceInMinutes <= minutes;
}


export const getUuidShort = (uuid: string): string => {
    return uuid.substring(0, 8)
}

export const getCamelCase = (str: string) => {
    return str.toLowerCase().replace(/[^a-zA-Z0-9]+(.)/g, (m, chr) => chr.toUpperCase());
}

export const sanitizeMongoDatabaseName = (name: string) => {
    
    // Explicit replacement according to MongoDB requirements:
    // https://www.mongodb.com/docs/manual/reference/limits/

    return getCamelCase(name.toLowerCase()
            .replaceAll(" ", "_")
            .replaceAll("/", "_")
            .replaceAll("\\", "_")
            .replaceAll(".", "_")
            .replaceAll('"', "_")
            .replaceAll('$', "_")
    )
}

export const formatAsPercentage = (num: number | string | null): string => {
    if (num === null) {
        return "-"
    }
    if (typeof num == "string") {
        num = +num
    }
    return new Intl.NumberFormat('default', {
        style: 'percent',
        minimumFractionDigits: 2,
        maximumFractionDigits: 2,
    }).format(num / 100);
}

export const formatAsThousands = (num: number | null): string => {
    if (num === null) {
        return "-"
    }
    return new Intl.NumberFormat('en-US').format(num);
}


export const baseTags = (tags: string[][] | undefined, unique: boolean = false, base: string[] = [FileTag.DNA, FileTag.RNA]): string[] => {
    if (tags === undefined) {
        return []
    }
    let baseTags = tags.map(t => {
        return base.map(b => {
            if (t.includes(b)) {
                return b
            }
            return null
        }).flatMap(f => f ? [f] : [])
    }).flat()

    if (unique){
        return baseTags.filter((x, i, a) => a.indexOf(x) == i)
    }

    return baseTags
}

export const baseTagsFromModel = (model: Cerebro, unique: boolean = false, base: string[] = [FileTag.DNA, FileTag.RNA]): string[] => {
    
    let baseTags = model.sample.tags.map(tag => {
        if (base.includes(tag)) {
            return tag
        } else {
            return null
        }
        }).flatMap(f => f ? [f] : []).flat()

    if (unique){
        return [...new Set(baseTags)]
    }
    return baseTags
}

export const getUserRoles = (adminRole: boolean, reportRole: boolean, dataRole: boolean, botRole: boolean): Role[] => {

    let roles = [Role.User];

    if (adminRole) {
        roles.push(Role.Admin)
    };
    if (reportRole) {
        roles.push(Role.Report)
    };
    if (dataRole) {
        roles.push(Role.Data)
    };
    if (botRole) {
        roles.push(Role.Bot)
    };

    return roles;
}

export const calculateOutput = (other_reads: number | null, phage_reads: number | null, qc_reads: number | null): number | null => {
    if (other_reads !== null && phage_reads !== null) {
        return Number(other_reads) - Number(phage_reads)
    } else if (other_reads !== null && phage_reads === null) {
        return other_reads
    } else if (other_reads === null && phage_reads !== null) {
        return phage_reads
    } else {
        return qc_reads
    }
}

export const calculateTotalQuality = (low_complexity_reads: number | null, qc_reads: number | null): number | null => {
    if (low_complexity_reads !== null && qc_reads !== null) {
        return Number(low_complexity_reads) + Number(qc_reads)
    } else if (low_complexity_reads === null && qc_reads !== null) {
        return qc_reads
    } else if (low_complexity_reads !== null && qc_reads === null) {
        return low_complexity_reads
    } else {
        return null
    }
}


   
export const getInitials = (name: string): string => {
    let nameSplit: string[] = name.split(" ")

    let initials: string = "";
    if (nameSplit.length > 1) {
        initials = `${nameSplit[0][0]}${nameSplit.slice(-1)[0]}`
    } else if (nameSplit.length == 1) {
        initials = `${name[0]}`
    } else {
        initials = `AU`
    }
    return initials;
}




export const ERCC_CONCENTRATIONS: Map<string, number> = new Map([
    ["ERCC-00130", 30000],
    ["ERCC-00004", 7500],
    ["ERCC-00136", 1875],
    ["ERCC-00108", 937.5],
    ["ERCC-00116", 468.75],
    ["ERCC-00092", 234.375],
    ["ERCC-00095", 117.1875],
    ["ERCC-00131", 117.1875],
    ["ERCC-00062", 58.59375],
    ["ERCC-00019", 29.296875],
    ["ERCC-00144", 29.296875],
    ["ERCC-00170", 14.6484375],
    ["ERCC-00154", 7.32421875],
    ["ERCC-00085", 7.32421875],
    ["ERCC-00028", 3.66210938],
    ["ERCC-00033", 1.83105469],
    ["ERCC-00134", 1.83105469],
    ["ERCC-00147", 0.91552734],
    ["ERCC-00097", 0.45776367],
    ["ERCC-00156", 0.45776367],
    ["ERCC-00123", 0.22888184],
    ["ERCC-00017", 0.11444092],
    ["ERCC-00083", 0.02861023],
    ["ERCC-00096", 15000],
    ["ERCC-00171", 3750],
    ["ERCC-00009", 937.5],
    ["ERCC-00042", 468.75],
    ["ERCC-00060", 234.375],
    ["ERCC-00035", 117.1875],
    ["ERCC-00025", 58.59375],
    ["ERCC-00051", 58.59375],
    ["ERCC-00053", 29.296875],
    ["ERCC-00148", 14.6484375],
    ["ERCC-00126", 14.6484375],
    ["ERCC-00034", 7.32421875],
    ["ERCC-00150", 3.66210938],
    ["ERCC-00067", 3.66210938],
    ["ERCC-00031", 1.83105469],
    ["ERCC-00109", 0.91552734],
    ["ERCC-00073", 0.91552734],
    ["ERCC-00158", 0.45776367],
    ["ERCC-00104", 0.22888184],
    ["ERCC-00142", 0.22888184],
    ["ERCC-00138", 0.11444092],
    ["ERCC-00117", 0.05722046],
    ["ERCC-00075", 0.01430512],
    ["ERCC-00074", 15000],
    ["ERCC-00113", 3750],
    ["ERCC-00145", 937.5],
    ["ERCC-00111", 468.75],
    ["ERCC-00076", 234.375],
    ["ERCC-00044", 117.1875],
    ["ERCC-00162", 58.59375],
    ["ERCC-00071", 58.59375],
    ["ERCC-00084", 29.296875],
    ["ERCC-00099", 14.6484375],
    ["ERCC-00054", 14.6484375],
    ["ERCC-00157", 7.32421875],
    ["ERCC-00143", 3.66210938],
    ["ERCC-00039", 3.66210938],
    ["ERCC-00058", 1.83105469],
    ["ERCC-00120", 0.91552734],
    ["ERCC-00040", 0.91552734],
    ["ERCC-00164", 0.45776367],
    ["ERCC-00024", 0.22888184],
    ["ERCC-00016", 0.22888184],
    ["ERCC-00012", 0.11444092],
    ["ERCC-00098", 0.05722046],
    ["ERCC-00057", 0.01430512],
    ["ERCC-00002", 15000],
    ["ERCC-00046", 3750],
    ["ERCC-00003", 937.5],
    ["ERCC-00043", 468.75],
    ["ERCC-00022", 234.375],
    ["ERCC-00112", 117.1875],
    ["ERCC-00165", 58.59375],
    ["ERCC-00079", 58.59375],
    ["ERCC-00078", 29.296875],
    ["ERCC-00163", 14.6484375],
    ["ERCC-00059", 14.6484375],
    ["ERCC-00160", 7.32421875],
    ["ERCC-00014", 3.66210938],
    ["ERCC-00077", 3.66210938],
    ["ERCC-00069", 1.83105469],
    ["ERCC-00137", 0.91552734],
    ["ERCC-00013", 0.91552734],
    ["ERCC-00168", 0.45776367],
    ["ERCC-00041", 0.22888184],
    ["ERCC-00081", 0.22888184],
    ["ERCC-00086", 0.11444092],
    ["ERCC-00061", 0.05722046],
    ["ERCC-00048", 0.01430512]
]);

export const ERCC_GROUPS: Map<string, string> = new Map([
    ["ERCC-00130", "A"],
    ["ERCC-00004", "A"],
    ["ERCC-00136", "A"],
    ["ERCC-00108", "A"],
    ["ERCC-00116", "A"],
    ["ERCC-00092", "A"],
    ["ERCC-00095", "A"],
    ["ERCC-00131", "A"],
    ["ERCC-00062", "A"],
    ["ERCC-00019", "A"],
    ["ERCC-00144", "A"],
    ["ERCC-00170", "A"],
    ["ERCC-00154", "A"],
    ["ERCC-00085", "A"],
    ["ERCC-00028", "A"],
    ["ERCC-00033", "A"],
    ["ERCC-00134", "A"],
    ["ERCC-00147", "A"],
    ["ERCC-00097", "A"],
    ["ERCC-00156", "A"],
    ["ERCC-00123", "A"],
    ["ERCC-00017", "A"],
    ["ERCC-00083", "A"],
    ["ERCC-00096", "B"],
    ["ERCC-00171", "B"],
    ["ERCC-00009", "B"],
    ["ERCC-00042", "B"],
    ["ERCC-00060", "B"],
    ["ERCC-00035", "B"],
    ["ERCC-00025", "B"],
    ["ERCC-00051", "B"],
    ["ERCC-00053", "B"],
    ["ERCC-00148", "B"],
    ["ERCC-00126", "B"],
    ["ERCC-00034", "B"],
    ["ERCC-00150", "B"],
    ["ERCC-00067", "B"],
    ["ERCC-00031", "B"],
    ["ERCC-00109", "B"],
    ["ERCC-00073", "B"],
    ["ERCC-00158", "B"],
    ["ERCC-00104", "B"],
    ["ERCC-00142", "B"],
    ["ERCC-00138", "B"],
    ["ERCC-00117", "B"],
    ["ERCC-00075", "B"],
    ["ERCC-00074", "C"],
    ["ERCC-00113", "C"],
    ["ERCC-00145", "C"],
    ["ERCC-00111", "C"],
    ["ERCC-00076", "C"],
    ["ERCC-00044", "C"],
    ["ERCC-00162", "C"],
    ["ERCC-00071", "C"],
    ["ERCC-00084", "C"],
    ["ERCC-00099", "C"],
    ["ERCC-00054", "C"],
    ["ERCC-00157", "C"],
    ["ERCC-00143", "C"],
    ["ERCC-00039", "C"],
    ["ERCC-00058", "C"],
    ["ERCC-00120", "C"],
    ["ERCC-00040", "C"],
    ["ERCC-00164", "C"],
    ["ERCC-00024", "C"],
    ["ERCC-00016", "C"],
    ["ERCC-00012", "C"],
    ["ERCC-00098", "C"],
    ["ERCC-00057", "C"],
    ["ERCC-00002", "D"],
    ["ERCC-00046", "D"],
    ["ERCC-00003", "D"],
    ["ERCC-00043", "D"],
    ["ERCC-00022", "D"],
    ["ERCC-00112", "D"],
    ["ERCC-00165", "D"],
    ["ERCC-00079", "D"],
    ["ERCC-00078", "D"],
    ["ERCC-00163", "D"],
    ["ERCC-00059", "D"],
    ["ERCC-00160", "D"],
    ["ERCC-00014", "D"],
    ["ERCC-00077", "D"],
    ["ERCC-00069", "D"],
    ["ERCC-00137", "D"],
    ["ERCC-00013", "D"],
    ["ERCC-00168", "D"],
    ["ERCC-00041", "D"],
    ["ERCC-00081", "D"],
    ["ERCC-00086", "D"],
    ["ERCC-00061", "D"],
    ["ERCC-00048", "D"]
]);

export const PATHOGENS: string[] = [
    'Abiotrophia defectiva',
    'Achromobacter xylosoxidans',
    'Acinetobacter baumannii',
    'Actinomadura madurae',
    'Actinomadura pelletieri',
    'Actinomyces gerencseriae',
    'Actinomyces israelii',
    'Actinotignum schaalii',
    'Aerococcus sanguinicola',
    'Aerococcus urinae',
    'Aeromonas caviae',
    'Aeromonas hydrophila',
    'Aeromonas veronii',
    'Aggregatibacter actinomycetemcomitans',
    'Aggregatibacter aphrophilus',
    'Agrobacterium tumefaciens',
    'Aliarcobacter butzleri',
    'Anaerococcus prevotii',
    'Anaplasma phagocytophilum',
    'Arcanobacterium haemolyticum',
    'Bacillus anthracis',
    'Bacillus cereus',
    'Bacillus licheniformis',
    'Bartonella bacilliformis',
    'Bartonella clarridgeiae',
    'Bartonella elizabethae',
    'Bartonella grahamii',
    'Bartonella henselae',
    'Bartonella quintana',
    'Bordetella bronchiseptica',
    'Bordetella parapertussis',
    'Bordetella pertussis',
    'Borrelia duttonii',
    'Borrelia hermsii',
    'Borrelia miyamotoi',
    'Borrelia parkeri',
    'Borrelia recurrentis',
    'Borrelia turicatae',
    'Borreliella burgdorferi',
    'Borreliella mayonii',
    'Brachyspira hyodysenteriae',
    'Brevibacillus brevis',
    'Brucella abortus',
    'Brucella canis',
    'Brucella melitensis',
    'Brucella suis',
    'Burkholderia cenocepacia',
    'Burkholderia cepacia',
    'Burkholderia mallei',
    'Burkholderia pseudomallei',
    'Campylobacter coli',
    'Campylobacter curvus',
    'Campylobacter fetus',
    'Campylobacter jejuni',
    'Campylobacter lari',
    'Cardiobacterium hominis',
    'Chlamydia pneumoniae',
    'Chlamydia psittaci',
    'Chlamydia trachomatis',
    'Citrobacter freundii',
    'Clostridioides difficile',
    'Clostridium botulinum',
    'Clostridium perfringens',
    'Clostridium septicum',
    'Clostridium tetani',
    'Corynebacterium diphtheriae',
    'Corynebacterium minutissimum',
    'Corynebacterium striatum',
    'Corynebacterium ulcerans',
    'Coxiella burnetii',
    'Cronobacter sakazakii',
    'Cronobacter turicensis',
    'Edwardsiella tarda',
    'Ehrlichia chaffeensis',
    'Ehrlichia ewingii',
    'Ehrlichia muris',
    'Eikenella corrodens',
    'Elizabethkingia anophelis',
    'Elizabethkingia meningoseptica',
    'Enterobacter cloacae',
    'Enterococcus avium',
    'Enterococcus faecalis',
    'Enterococcus faecium',
    'Enterococcus gallinarum',
    'Erysipelothrix rhusiopathiae',
    'Escherichia coli',
    'Finegoldia magna',
    'Francisella tularensis',
    'Fusobacterium necrophorum',
    'Fusobacterium nucleatum',
    'Gardnerella vaginalis',
    'Haemophilus aegyptius',
    'Haemophilus influenzae',
    'Haemophilus parainfluenzae',
    'Helicobacter pylori',
    'Kingella denitrificans',
    'Kingella kingae',
    'Klebsiella aerogenes',
    'Klebsiella granulomatis',
    'Klebsiella oxytoca',
    'Klebsiella pneumoniae',
    'Legionella bozemanae',
    'Legionella micdadei',
    'Legionella pneumophila',
    'Leptospira alexanderi',
    'Leptospira alstonii',
    'Leptospira borgpetersenii',
    'Leptospira interrogans',
    'Leptospira kirschneri',
    'Leptospira kmetyi',
    'Leptospira mayottensis',
    'Leptospira noguchii',
    'Leptospira santarosai',
    'Leptospira weilii',
    'Listeria monocytogenes',
    'Malacoplasma penetrans',
    'Metamycoplasma hominis',
    'Metamycoplasma salivarium',
    'Methylorubrum extorquens',
    'Moraxella catarrhalis',
    'Moraxella lacunata',
    'Moraxella nonliquefaciens',
    'Morganella morganii',
    'Mycobacterium avium',
    'Mycobacterium intracellulare',
    'Mycobacterium kansasii',
    'Mycobacterium leprae',
    'Mycobacterium lepromatosis',
    'Mycobacterium malmoense',
    'Mycobacterium marinum',
    'Mycobacterium scrofulaceum',
    'Mycobacterium simiae',
    'Mycobacterium szulgai',
    'Mycobacterium tuberculosis',
    'Mycobacterium ulcerans',
    'Mycobacterium xenopi',
    'Mycobacteroides abscessus',
    'Mycobacteroides chelonae',
    'Mycolicibacter terrae',
    'Mycolicibacterium fortuitum',
    'Mycolicibacterium smegmatis',
    'Mycolicibacterium wolinskyi',
    'Mycoplasmoides genitalium',
    'Mycoplasmoides pneumoniae',
    'Mycoplasmopsis caviae',
    'Mycoplasmopsis fermentans',
    'Neisseria gonorrhoeae',
    'Neisseria meningitidis',
    'Neorickettsia sennetsu',
    'Nocardia asteroides',
    'Nocardia brasiliensis',
    'Nocardia farcinica',
    'Nocardia nova',
    'Nocardia otitidiscaviarum',
    'Orientia tsutsugamushi',
    'Paraclostridium bifermentans',
    'Parvimonas micra',
    'Pasteurella multocida',
    'Peptoniphilus asaccharolyticus',
    'Peptostreptococcus anaerobius',
    'Photobacterium damselae',
    'Plesiomonas shigelloides',
    'Porphyromonas gingivalis',
    'Prescottella equi',
    'Prevotella bivia',
    'Prevotella intermedia',
    'Prevotella melaninogenica',
    'Proteus mirabilis',
    'Proteus penneri',
    'Proteus vulgaris',
    'Providencia alcalifaciens',
    'Providencia rettgeri',
    'Pseudomonas aeruginosa',
    'Pseudomonas fluorescens',
    'Pseudomonas putida',
    'Ralstonia pickettii',
    'Rickettsia africae',
    'Rickettsia akari',
    'Rickettsia australis',
    'Rickettsia canadensis',
    'Rickettsia conorii',
    'Rickettsia felis',
    'Rickettsia helvetica',
    'Rickettsia honei',
    'Rickettsia japonica',
    'Rickettsia montanensis',
    'Rickettsia parkeri',
    'Rickettsia prowazekii',
    'Rickettsia rickettsii',
    'Rickettsia sibirica',
    'Rickettsia typhi',
    'Rothia dentocariosa',
    'Salmonella bongori',
    'Salmonella enterica',
    'Serratia marcescens',
    'Shewanella algae',
    'Shigella boydii',
    'Shigella dysenteriae',
    'Shigella flexneri',
    'Shigella sonnei',
    'Staphylococcus aureus',
    'Staphylococcus capitis',
    'Staphylococcus epidermidis',
    'Staphylococcus haemolyticus',
    'Staphylococcus lugdunensis',
    'Staphylococcus pseudintermedius',
    'Staphylococcus saccharolyticus',
    'Staphylococcus saprophyticus',
    'Stenotrophomonas maltophilia',
    'Streptobacillus moniliformis',
    'Streptococcus agalactiae',
    'Streptococcus anginosus',
    'Streptococcus criceti',
    'Streptococcus dysgalactiae',
    'Streptococcus equi',
    'Streptococcus equinus',
    'Streptococcus gordonii',
    'Streptococcus iniae',
    'Streptococcus mitis',
    'Streptococcus mutans',
    'Streptococcus oralis',
    'Streptococcus pneumoniae',
    'Streptococcus pyogenes',
    'Streptococcus sanguinis',
    'Streptococcus sobrinus',
    'Streptococcus suis',
    'Suttonella indologenes',
    'Tannerella forsythia',
    'Treponema denticola',
    'Treponema pallidum',
    'Treponema vincentii',
    'Tropheryma whipplei',
    'Trueperella pyogenes',
    'Ureaplasma parvum',
    'Ureaplasma urealyticum',
    'Vibrio alginolyticus',
    'Vibrio cholerae',
    'Vibrio parahaemolyticus',
    'Vibrio vulnificus',
    'Yersinia enterocolitica',
    'Yersinia pestis',
    'Yersinia pseudotuberculosis',
    '[Haemophilus] ducreyi',
    'Acanthamoeba astronyxis',
    'Acanthamoeba castellanii',
    'Acanthamoeba culbertsoni',
    'Acanthamoeba divionensis',
    'Acanthamoeba hatchetti',
    'Acanthamoeba lenticulata',
    'Acanthamoeba lugdunensis',
    'Acanthamoeba polyphaga',
    'Acanthamoeba rhysodes',
    'Adenocephalus pacificus',
    'Alternaria alternata',
    'Alternaria brassicicola',
    'Ancylostoma braziliense',
    'Ancylostoma caninum',
    'Ancylostoma ceylanicum',
    'Ancylostoma duodenale',
    'Angiostrongylus cantonensis',
    'Angiostrongylus costaricensis',
    'Anisakis simplex',
    'Anncaliia algerae',
    'Apophysomyces variabilis',
    'Ascaris lumbricoides',
    'Ascaris suum',
    'Aspergillus flavus',
    'Aspergillus fumigatus',
    'Aspergillus nidulans',
    'Aspergillus niger',
    'Aspergillus terreus',
    'Babesia divergens',
    'Babesia duncani',
    'Babesia microti',
    'Babesia sp. venatorum',
    'Balamuthia mandrillaris',
    'Basidiobolus ranarum',
    'Baylisascaris procyonis',
    'Bertiella mucronata',
    'Bertiella studeri',
    'Blastocystis hominis',
    'Blastomyces dermatitidis',
    'Blastomyces gilchristii',
    'Blastomyces helicus',
    'Brugia malayi',
    'Brugia pahangi',
    'Brugia timori',
    'Candida albicans',
    'Candida dubliniensis',
    'Candida parapsilosis',
    'Candida tropicalis',
    'Capillaria hepatica',
    'Chrysomya bezziana',
    'Cladophialophora bantiana',
    'Cladophialophora boppii',
    'Cladophialophora carrionii',
    'Clavispora lusitaniae',
    'Clinostomum complanatum',
    'Clonorchis sinensis',
    'Coccidioides immitis',
    'Coccidioides posadasii',
    'Cochliomyia hominivorax',
    'Colletotrichum coccodes',
    'Colletotrichum gloeosporioides',
    'Conidiobolus coronatus',
    'Conidiobolus incongruus',
    'Contracaecum osculatum',
    'Cryptococcus gattii',
    'Cryptococcus neoformans',
    'Cryptosporidium hominis',
    'Cryptosporidium parvum',
    'Curvularia australiensis',
    'Curvularia hawaiiensis',
    'Curvularia lunata',
    'Curvularia papendorfii',
    'Curvularia spicifera',
    'Cyclospora cayetanensis',
    'Cystoisospora belli',
    'Dermanyssus gallinae',
    'Dermatobia hominis',
    'Dicrocoelium dendriticum',
    'Dientamoeba fragilis',
    'Dioctophyme renale',
    'Diphyllobothrium latum',
    'Dirofilaria immitis',
    'Dirofilaria repens',
    'Dracunculus medinensis',
    'Echinococcus granulosus',
    'Echinococcus multilocularis',
    'Echinococcus oligarthrus',
    'Echinococcus vogeli',
    'Echinolaelaps echidninus',
    'Emmonsia crescens',
    'Encephalitozoon cuniculi',
    'Encephalitozoon hellem',
    'Encephalitozoon intestinalis',
    'Entamoeba histolytica',
    'Enterobius vermicularis',
    'Enterocytozoon bieneusi',
    'Epidermophyton floccosum',
    'Exophiala dermatitidis',
    'Fasciola gigantica',
    'Fasciola hepatica',
    'Fasciolopsis buski',
    'Fonsecaea compacta',
    'Fonsecaea pedrosoi',
    'Fusarium fujikuroi',
    'Fusarium oxysporum',
    'Fusarium solani',
    'Giardia intestinalis',
    'Gnathostoma hispidum',
    'Gnathostoma spinigerum',
    'Halicephalobus gingivalis',
    'Heterophyes heterophyes',
    'Histoplasma capsulatum',
    'Hymenolepis diminuta',
    'Leishmania aethiopica',
    'Leishmania amazonensis',
    'Leishmania braziliensis',
    'Leishmania chagasi',
    'Leishmania colombiensis',
    'Leishmania donovani',
    'Leishmania guyanensis',
    'Leishmania infantum',
    'Leishmania major',
    'Leishmania martiniquensis',
    'Leishmania mexicana',
    'Leishmania naiffi',
    'Leishmania panamensis',
    'Leishmania peruviana',
    'Leishmania tropica',
    'Lichtheimia corymbifera',
    'Lichtheimia ramosa',
    'Linguatula serrata',
    'Loa loa',
    'Lomentospora prolificans',
    'Madurella mycetomatis',
    'Mansonella ozzardi',
    'Mansonella perstans',
    'Mansonella streptocerca',
    'Metagonimus yokogawai',
    'Moniliformis moniliformis',
    'Mucor circinelloides',
    'Naegleria fowleri',
    'Nakaseomyces glabratus',
    'Necator americanus',
    'Neobalantidium coli',
    'Onchocerca volvulus',
    'Opisthorchis felineus',
    'Opisthorchis sinensis',
    'Opisthorchis viverrini',
    'Ornithonyssus bacoti',
    'Ornithonyssus bursa',
    'Ornithonyssus sylviarum',
    'Paracapillaria philippinensis',
    'Paracoccidioides brasiliensis',
    'Paragonimus africanus',
    'Paragonimus heterotremus',
    'Paragonimus kellicotti',
    'Paragonimus mexicanus',
    'Paragonimus siamensis',
    'Paragonimus skrjabini',
    'Paragonimus uterobilateralis',
    'Paragonimus westermani',
    'Pediculus humanus',
    'Phialophora verrucosa',
    'Pichia kudriavzevii',
    'Plasmodium falciparum',
    'Plasmodium knowlesi',
    'Plasmodium malariae',
    'Plasmodium ovale',
    'Plasmodium vivax',
    'Pneumocystis jirovecii',
    'Pseudoterranova decipiens',
    'Pthirus pubis',
    'Pythium insidiosum',
    'Rhinocladiella mackenziei',
    'Rhinosporidium seeberi',
    'Rhizomucor pusillus',
    'Rhizopus arrhizus',
    'Rhizopus microsporus',
    'Rodentolepis nana',
    'Saccharomyces cerevisiae',
    'Saksenaea erythrospora',
    'Saksenaea vasiformis',
    'Sappinia diploidea',
    'Sappinia pedata',
    'Sarcocystis hominis',
    'Sarcocystis suihominis',
    'Sarcoptes scabiei',
    'Scedosporium apiospermum',
    'Scedosporium boydii',
    'Schistosoma guineensis',
    'Schistosoma haematobium',
    'Schistosoma intercalatum',
    'Schistosoma japonicum',
    'Schistosoma malayensis',
    'Schistosoma mansoni',
    'Schistosoma mekongi',
    'Scopulariopsis brevicaulis',
    'Spirometra erinaceieuropaei',
    'Sporothrix schenckii',
    'Strongyloides stercoralis',
    'Taenia multiceps',
    'Taenia saginata',
    'Taenia solium',
    'Talaromyces marneffei',
    'Thelazia callipaeda',
    'Toxocara canis',
    'Toxocara cati',
    'Toxoplasma gondii',
    'Trachipleistophora hominis',
    'Trematosphaeria grisea',
    'Trichinella britovi',
    'Trichinella nativa',
    'Trichinella nelsoni',
    'Trichinella pseudospiralis',
    'Trichinella spiralis',
    'Trichomonas vaginalis',
    'Trichophyton rubrum',
    'Trichuris trichiura',
    'Trichuris vulpis',
    'Trypanosoma brucei',
    'Trypanosoma cruzi',
    'Tubulinosema acridophagus',
    'Tunga penetrans',
    'Uncinaria stenocephala',
    'Vittaforma corneae',
    'Wuchereria bancrofti',
    '[Candida] auris',
    'Aichivirus A',
    'Akhmeta virus',
    'Alphainfluenzavirus influenzae',
    'Alphainfluenzavirus influenzae',
    'Alphapapillomavirus 4',
    'Alphapapillomavirus 7',
    'Alphapapillomavirus 9',
    'Alphapolyomavirus octihominis',
    'Alphapolyomavirus octihominis',
    'Alphapolyomavirus quintihominis',
    'Alphapolyomavirus quintihominis',
    'Bandavirus bhanjanagarense',
    'Bandavirus bhanjanagarense',
    'Bandavirus dabieense',
    'Banna virus',
    'Barmah Forest virus',
    'Bayou orthohantavirus',
    'Betacoronavirus 1',
    'Betainfluenzavirus influenzae',
    'Betainfluenzavirus influenzae',
    'Betapolyomavirus hominis',
    'Betapolyomavirus quartihominis',
    'Betapolyomavirus secuhominis',
    'Bhanja bandavirus',
    'Black Creek Canal orthohantavirus',
    'Borna disease virus',
    'Bufavirus-3',
    'Bunyamwera orthobunyavirus',
    'Bwamba orthobunyavirus',
    'Cache Valley orthobunyavirus',
    'Cali mammarenavirus',
    'California encephalitis orthobunyavirus',
    'Cardiovirus A',
    'Chapare mammarenavirus',
    'Chikungunya virus',
    'Colorado tick fever coltivirus',
    'Colorado tick fever coltivirus',
    'Cowpox virus',
    'Coxsackievirus',
    'Crimean-Congo hemorrhagic fever orthonairovirus',
    'Cytomegalovirus humanbeta5',
    'Cytomegalovirus humanbeta5',
    'Deltapolyomavirus decihominis',
    'Deltapolyomavirus decihominis',
    'Dobrava-Belgrade orthohantavirus',
    'Dugbe orthonairovirus',
    'Eastern equine encephalitis virus',
    'Echovirus',
    'Enterovirus A',
    'Enterovirus B',
    'Enterovirus C',
    'Enterovirus D',
    'Erythroparvovirus primate1',
    'European bat lyssavirus',
    'Everglades virus',
    'Gammainfluenzavirus influenzae',
    'Gammainfluenzavirus influenzae',
    'Guanarito mammarenavirus',
    'Guaroa orthobunyavirus',
    'Hantaan orthohantavirus',
    'Heartland bandavirus',
    'Hendra henipavirus',
    'Hepacivirus C',
    'Hepatitis B virus',
    'Hepatitis delta virus',
    'Hepatovirus A',
    'Huaiyangshan banyangvirus',
    'Human adenovirus sp.',
    'Human astrovirus',
    'Human betaherpesvirus 6',
    'Human bocavirus',
    'Human coronavirus 229E',
    'Human coronavirus HKU1',
    'Human coronavirus NL63',
    'Human immunodeficiency virus',
    'Human mastadenovirus A',
    'Human mastadenovirus B',
    'Human mastadenovirus C',
    'Human mastadenovirus D',
    'Human mastadenovirus F',
    'Human parechovirus',
    'Human rhinovirus sp.',
    'Human torovirus',
    'Ilesha orthobunyavirus',
    'Jamestown Canyon orthobunyavirus',
    'La Crosse orthobunyavirus',
    'Lassa mammarenavirus',
    'Lujo mammarenavirus',
    'Lymphocryptovirus humangamma4',
    'Lymphocryptovirus humangamma4',
    'Lymphocytic choriomeningitis mammarenavirus',
    'Lyssavirus australis',
    'Lyssavirus duvenhage',
    'Lyssavirus mokola',
    'Lyssavirus rabies',
    'Machupo mammarenavirus',
    'Mammalian orthoreovirus',
    'Mammalian rubulavirus 5',
    'Mammarenavirus brazilense',
    'Mammarenavirus juninense',
    'Marburg marburgvirus',
    'Mayaro virus',
    'Metapneumovirus hominis',
    'Metapneumovirus hominis',
    'Middelburg virus',
    'Middle East respiratory syndrome-related coronavirus',
    'Molluscum contagiosum virus',
    'Monkeypox virus',
    'Morbillivirus hominis',
    'Morbillivirus hominis',
    'Mucambo virus',
    'Mumps rubulavirus',
    'Mupapillomavirus 1',
    'Nairobi sheep disease orthonairovirus',
    'Naples phlebovirus',
    'Nipah henipavirus',
    'Norwalk virus',
    'Onyong-nyong virus',
    'Orf virus',
    'Oropouche orthobunyavirus',
    'Orthoavulavirus javaense',
    'Orthoavulavirus javaense',
    'Orthobunyavirus bunyamweraense',
    'Orthobunyavirus cacheense',
    'Orthobunyavirus guaroaense',
    'Orthobunyavirus ileshaense',
    'Orthobunyavirus khatangaense',
    'Orthobunyavirus lacrosseense',
    'Orthobunyavirus oropoucheense',
    'Orthoebolavirus bundibugyoense',
    'Orthoebolavirus bundibugyoense',
    'Orthoebolavirus sudanense',
    'Orthoebolavirus sudanense',
    'Orthoebolavirus zairense',
    'Orthoebolavirus zairense',
    'Orthoflavivirus denguei',
    'Orthoflavivirus denguei',
    'Orthoflavivirus encephalitidis',
    'Orthoflavivirus encephalitidis',
    'Orthoflavivirus flavi',
    'Orthoflavivirus flavi',
    'Orthoflavivirus ilheusense',
    'Orthoflavivirus ilheusense',
    'Orthoflavivirus japonicum',
    'Orthoflavivirus japonicum',
    'Orthoflavivirus kyasanurense',
    'Orthoflavivirus kyasanurense',
    'Orthoflavivirus louisense',
    'Orthoflavivirus louisense',
    'Orthoflavivirus loupingi',
    'Orthoflavivirus loupingi',
    'Orthoflavivirus murrayense',
    'Orthoflavivirus murrayense',
    'Orthoflavivirus nilense',
    'Orthoflavivirus nilense',
    'Orthoflavivirus omskense',
    'Orthoflavivirus omskense',
    'Orthoflavivirus powassanense',
    'Orthoflavivirus powassanense',
    'Orthoflavivirus wesselsbronense',
    'Orthoflavivirus wesselsbronense',
    'Orthoflavivirus zikaense',
    'Orthoflavivirus zikaense',
    'Orthohantavirus andesense',
    'Orthohantavirus sinnombreense',
    'Orthohantavirus sinnombreense',
    'Orthonairovirus nairobiense',
    'Orthopneumovirus hominis',
    'Orthopneumovirus hominis',
    'Orthorubulavirus hominis',
    'Orthorubulavirus hominis',
    'Orthorubulavirus laryngotracheitidis',
    'Orthorubulavirus laryngotracheitidis',
    'Orungo virus',
    'Parechovirus A',
    'Parechovirus B',
    'Paslahepevirus balayani',
    'Pegivirus C',
    'Pegivirus hominis',
    'Phlebovirus napoliense',
    'Phlebovirus riftense',
    'Phlebovirus siciliaense',
    'Phlebovirus siciliaense',
    'Phlebovirus toscanaense',
    'Primate T-lymphotropic virus 1',
    'Primate T-lymphotropic virus 2',
    'Primate erythroparvovirus 1',
    'Pseudocowpox virus',
    'Punta Toro phlebovirus',
    'Puumala orthohantavirus',
    'Respirovirus laryngotracheitidis',
    'Respirovirus laryngotracheitidis',
    'Respirovirus pneumoniae',
    'Respirovirus pneumoniae',
    'Rhadinovirus humangamma8',
    'Rhadinovirus humangamma8',
    'Rhinovirus A',
    'Rhinovirus C',
    'Rift Valley fever phlebovirus',
    'Rosavirus A',
    'Roseolovirus humanbeta6a',
    'Roseolovirus humanbeta6a',
    'Roseolovirus humanbeta6b',
    'Roseolovirus humanbeta6b',
    'Roseolovirus humanbeta7',
    'Roseolovirus humanbeta7',
    'Ross River virus',
    'Rotavirus A',
    'Rotavirus B',
    'Rotavirus C',
    'Rubivirus rubellae',
    'Rubivirus rubellae',
    'Russian Spring-Summer encephalitis virus',
    'Salivirus A',
    'Sandfly fever Naples phlebovirus',
    'Sapporo virus',
    'Semliki Forest virus',
    'Seoul orthohantavirus',
    'Severe acute respiratory syndrome-related coronavirus',
    'Sicilian phlebovirus',
    'Simplexvirus humanalpha1',
    'Simplexvirus humanalpha1',
    'Simplexvirus humanalpha2',
    'Simplexvirus humanalpha2',
    'Simplexvirus macacinealpha1',
    'Simplexvirus macacinealpha1',
    'Sin Nombre orthohantavirus',
    'Sindbis virus',
    'Snowshoe hare orthobunyavirus',
    'Spondweni virus',
    'Tanapox virus',
    'Thogotovirus bourbonense',
    'Thogotovirus bourbonense',
    'Thogotovirus dhoriense',
    'Tonate virus',
    'Toscana phlebovirus',
    'Uukuniemi phlebovirus',
    'Uukuniemi uukuvirus',
    'Uukuvirus uukuniemiense',
    'Vaccinia virus',
    'Varicellovirus humanalpha3',
    'Varicellovirus humanalpha3',
    'Variola virus',
    'Venezuelan equine encephalitis virus',
    'Vesicular stomatitis virus',
    'Vesiculovirus chandipura',
    'Vesiculovirus isfahan',
    'Western equine encephalitis virus',
    'Whitewater Arroyo mammarenavirus',
    'Yaba monkey tumor virus',
    'unidentified human coronavirus'
]