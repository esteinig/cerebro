
import { FileTag, Role } from "./types";


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
