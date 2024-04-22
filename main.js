import 'bootstrap/dist/css/bootstrap.min.css';

const vcfForm = document.querySelector("#vcfForm");
const vcfInput = document.querySelector("#vcfInput");
const namesInput = document.querySelector("#namesInput");
const progressBar = document.querySelector("#progressBar");
const startButton = document.querySelector("#startButton");
const errorText = document.querySelector("#errorText");

vcfForm.addEventListener('submit', e => {
    e.preventDefault();
    startButton.innerText = 'Cancel';
    errorText.innerText = '';

    const sampleNames = namesInput.value.split(/\s/);

    const [vcfFile] = vcfInput.files;
    const reader = vcfFile.stream().getReader();

    startButton.addEventListener('click', async e => {
        e.preventDefault();
        await reader.cancel();
    }, { once: true });

    convert(reader, vcfFile.size, sampleNames)
        .then(blob => {
            downloadBlob(blob, "output.fasta");
        })
        .catch(err => {
            console.error(err);
            errorText.innerText = err.message;
        })
        .finally(() => {
            startButton.innerText = 'Start';
        });
});

function getReplacement(location, reference, alternates) {
    if (location === ".") return null;
    return location === "0" ? reference : alternates[location - 1];
}

async function convert(stream, totalSize, sampleNames) {
    function handleProgress(percentage) {
        return new Promise(res => {
            requestAnimationFrame(() => {
                progressBar.style.width = `${percentage * 100}%`;
                res();
            });
        });
    }

    const sampleCount = sampleNames.length;
    const results = [];

    sampleNames.push("ref");
    for (let si = 0; si < sampleNames.length; si++)
        results.push([]);

    let window = [];
    let windowEnd = 0;
    let headers = null;
    let prevChromosome = "";

    for await (const line of makeTextFileLineIterator(stream, totalSize, handleProgress)) {
        // Skip comments
        if (line.startsWith("##")) continue;

        const columns = line.split('\t');
        if (headers == null) {
            // Parse header
            const lastColumnIndex = columns.length - 1;
            const actualSampleCount = lastColumnIndex - columns.indexOf("FORMAT");
            if (sampleCount !== actualSampleCount)
                throw new Error(`Found ${actualSampleCount} samples, but ${sampleCount} names were given.`);

            headers = columns;
            continue;
        }

        const namedColumns = {};
        for (let i = 0; i < columns.length; i++)
            namedColumns[headers[i]] = columns[i];

        const chromosome = namedColumns["#CHROM"];
        if (chromosome !== prevChromosome) {
            console.log("Entering chromosome", chromosome);

            // Reset position state
            windowEnd = 0;
            prevChromosome = chromosome;
        }

        const position = parseInt(namedColumns["POS"]);
        if (position >= windowEnd) {
            for (const replacements of window) {
                let maxReplacementLength = 0;
                for (let si = 0; si < replacements.length; si++) {
                    const {length} = replacements[si] ??= "N";
                    if (length > maxReplacementLength)
                        maxReplacementLength = length;
                }

                for (let si = 0; si < replacements.length; si++) {
                    const replacement = replacements[si].padEnd(maxReplacementLength, "-");
                    results[si].push(replacement);
                }
            }

            window = [];
            windowEnd = position;
        }

        // Expand window to fit max referenced position
        const reference = namedColumns["REF"];
        const minWindowEnd = position + reference.length;
        for (let pos = windowEnd; pos < minWindowEnd; pos++)
            window.push(new Array(sampleNames.length));

        // Update window end position
        if (minWindowEnd > windowEnd)
            windowEnd = minWindowEnd;

        // Calculate position in window
        const windowStart = windowEnd - window.length;
        const windowPos = position - windowStart;

        const alternates = namedColumns["ALT"].split(',');
        const samples = columns.slice(-sampleCount);
        for (let si = 0; si < samples.length; si++) {
            if (window[windowPos][si] != null) continue;

            const sample = samples[si];
            const location = /^(\d+|\.)/.exec(sample)[0];
            const replacement = getReplacement(location, reference, alternates);
            if (replacement == null) continue;

            if (replacement.length >= reference.length) {
                window[windowPos][si] = replacement;
                continue;
            }

            for (let i = 0; i < replacement.length; i++)
                window[windowPos + i][si] = replacement[i];

            for (let i = replacement.length; i < reference.length; i++)
                window[windowPos + i][si] = '-';
        }

        for (let i = 0; i < reference.length; i++)
            window[windowPos + i][sampleCount] = reference[i];
    }

    for (let si = 0; si < sampleNames.length; si++)
        results[si] = `>${sampleNames[si]}\n${results[si].join("")}\n`;

    return new Blob(results, { type: 'text/plain' });
}

function downloadBlob(blob, filename){
    const a = document.createElement('a') // Create "a" element
    const url = URL.createObjectURL(blob) // Create an object URL from blob
    a.setAttribute('href', url) // Set "a" element link
    a.setAttribute('download', filename) // Set download filename
    a.click() // Start downloading
}

async function* asyncRead(reader) {
    while (true) {
        const { done, value } = await reader.read();
        yield value;
        if (done) break;
    }
}

async function* makeTextFileLineIterator(reader, totalSize, onProgress) {
    const utf8Decoder = new TextDecoder("utf-8");
    const decodeChunk = chunk =>
        chunk ? utf8Decoder.decode(chunk, { stream: true }) : "";

    let readSize = 0;
    let chunkText = "";
    for await (const chunk of asyncRead(reader)) {
        chunkText += decodeChunk(chunk);

        const re = /\r\n|\n|\r/gm;
        let result, startIndex = 0;
        while (result = re.exec(chunkText)) {
            yield chunkText.substring(startIndex, result.index);
            startIndex = re.lastIndex;
        }
        chunkText = chunkText.substring(startIndex);

        const chunkLength = chunk?.length ?? 0;
        readSize += chunkLength;

        await onProgress(readSize / totalSize);
    }

    if (readSize < totalSize)
        throw new Error("Operation was canceled");

    if (chunkText)
        yield chunkText;
}
