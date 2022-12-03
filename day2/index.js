const readInput = require('../libs/readInput')
const assert = require("assert")

// A = rock
// B = paper
// C = scissors
// X = Rock
// Y = Paper
// Z = Scissors

const [lost, draw, win] = [0, 3, 6]
const [rock, paper, scissors] = [1, 2, 3]


const map = {
    'A X': draw + rock, // r/r - draw rock
    'A Y': win + paper, // r/p - win paper
    'A Z': lost + scissors, // r/s - lost scissors

    'B X': lost + rock, // p/x - lost rock
    'B Y': draw + paper, // p/p - draw paper
    'B Z': win + scissors, // p/s - win scissors

    'C X': win + rock, // s/x - win rock
    'C Y': lost + paper, // s/p - lost paper
    'C Z': draw + scissors, // s/s - draw scissors
}

const testInput = [
    'A Y',
    'B X',
    'C Z',
]

const desiredMap = {
    Y: {A: 'X', B: 'Y', C: 'Z'}, // draw
    Z: {A: 'Y', B: 'Z', C: 'X'}, // win
    X: {A: 'Z', B: 'X', C: 'Y'} // lose
}

const basicSum = (input) => input.reduce((acc, line) => acc + map[line], 0)

const desiredSum = (input) => input.reduce((acc, line) => {
    const [opponent, result] = line.split(' ')
    return acc + map[`${opponent} ${desiredMap[result][opponent]}`]
}, 0)

const main = async () => {
    assert(basicSum(testInput) === 15, `basicSum failed, expected 15 but got ${basicSum(testInput)}`)
    assert(desiredSum(testInput) === 12, `desiredSum failed, expected 12 but got ${desiredSum(testInput)}`)
    const lines = await readInput(__dirname, 'input.txt')
    console.log(`Part 1 Result: ${basicSum(lines)}`)
    console.log(`Part 2 Result: ${desiredSum(lines)}`)
}

main().catch(console.error)
