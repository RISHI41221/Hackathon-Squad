"use strict";

const fs = require("fs");
const path = require("path");
const { execFileSync } = require("child_process");
const { parseInput, validateSolution } = require("./solver");

const SOLVER_PATH = path.join(__dirname, "hackathon_squad_solver.exe");
const SOLVER_TIMEOUT_MS = Number(process.env.HS_HARNESS_TIMEOUT_MS || 300000);

function formatRatio(ours, expected) {
  if (expected === 0) {
    return "n/a";
  }
  return `${((ours / expected) * 100).toFixed(2)}%`;
}

function parseSolverOutput(outputText, vertexCount) {
  const normalized = outputText.replace(/\r/g, "").trimEnd();
  const lines = normalized.length === 0 ? [] : normalized.split("\n");
  if (lines.length === 0) {
    throw new Error("Solver produced empty output");
  }

  const totalWeight = Number(lines[0].trim());
  if (!Number.isFinite(totalWeight)) {
    throw new Error(`Invalid total weight line: ${JSON.stringify(lines[0])}`);
  }

  const chosen = new Uint8Array(vertexCount);
  const pickedLine = lines.length >= 2 ? lines[1].trim() : "";
  const picked = pickedLine === ""
    ? []
    : pickedLine.split(/\s+/).map((token) => {
      const value = Number(token);
      if (!Number.isInteger(value)) {
        throw new Error(`Invalid picked vertex: ${JSON.stringify(token)}`);
      }
      return value;
    });

  for (const vertex of picked) {
    if (vertex < 1 || vertex > vertexCount) {
      throw new Error(`Picked vertex out of range: ${vertex}`);
    }
    chosen[vertex - 1] = 1;
  }

  return { totalWeight, picked, chosen };
}

function solveWithBinary(inputText, instance) {
  const stdout = execFileSync(SOLVER_PATH, [], {
    cwd: __dirname,
    input: inputText,
    encoding: "utf8",
    timeout: SOLVER_TIMEOUT_MS,
    maxBuffer: 16 * 1024 * 1024,
    windowsHide: true
  });

  return parseSolverOutput(stdout, instance.n);
}

function main() {
  const suiteDir = process.argv[2] || "C:\\Users\\RISHI\\Downloads\\test_suite\\test_suite";
  const files = fs.readdirSync(suiteDir)
    .filter((name) => /^input.*\.txt$/i.test(name))
    .sort((a, b) => a.localeCompare(b));

  let passCount = 0;
  let totalCases = 0;

  for (const name of files) {
    totalCases += 1;
    const inputPath = path.join(suiteDir, name);
    const expectedPath = path.join(suiteDir, name.replace(/^input/i, "output"));
    const inputText = fs.readFileSync(inputPath, "utf8");
    const instance = parseInput(inputText);
    const expectedWeight = Number(fs.readFileSync(expectedPath, "utf8").trim().split(/\s+/)[0]);

    try {
      const started = Date.now();
      const result = solveWithBinary(inputText, instance);
      const elapsedMs = Date.now() - started;
      const verifiedWeight = validateSolution(instance, result.chosen);
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
    } catch (error) {
      const reason = error && error.code === "ETIMEDOUT" ? "CRASH/TIMEOUT" : "CRASH/TIMEOUT";
      console.log([name, reason].join("\t"));
    }
  }

  console.log(`summary\tvalid_cases=${passCount}/${totalCases}`);
}

main();
