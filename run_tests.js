"use strict";

const fs = require("fs");
const path = require("path");
const { parseInput, solveInstance, validateSolution } = require("./solver");

function formatRatio(ours, expected) {
  if (expected === 0) {
    return "n/a";
  }
  return `${((ours / expected) * 100).toFixed(2)}%`;
}

function main() {
  const suiteDir = process.argv[2] || "C:\\Users\\RISHI\\Downloads\\test_suite\\test_suite";
  const files = fs.readdirSync(suiteDir)
    .filter((name) => /^input_.*\.txt$/i.test(name))
    .sort((a, b) => a.localeCompare(b));

  let passCount = 0;
  let totalCases = 0;

  for (const name of files) {
    totalCases += 1;
    const inputPath = path.join(suiteDir, name);
    const expectedPath = path.join(suiteDir, name.replace(/^input/i, "output"));
    const inputText = fs.readFileSync(inputPath, "utf8");
    const instance = parseInput(inputText);

    const started = Date.now();
    const result = solveInstance(instance);
    const elapsedMs = Date.now() - started;
    const verifiedWeight = validateSolution(instance, result.chosen);
    const expectedWeight = Number(fs.readFileSync(expectedPath, "utf8").trim().split(/\s+/)[0]);
    const ok = Math.round(verifiedWeight) === Math.round(result.totalWeight);
    if (ok) {
      passCount += 1;
    }

    console.log(
      [
        name,
        `valid=${ok ? "yes" : "no"}`,
        `ours=${Math.round(result.totalWeight)}`,
        `expected=${expectedWeight}`,
        `ratio=${formatRatio(result.totalWeight, expectedWeight)}`,
        `time_ms=${elapsedMs}`
      ].join("\t")
    );
  }

  console.log(`summary\tvalid_cases=${passCount}/${totalCases}`);
}

main();
