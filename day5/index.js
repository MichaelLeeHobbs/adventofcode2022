const readInput = require('../libs/readInput')
const assert = require("assert")
const fs = require("fs/promises")
const path = require("path")

const testInput = [
    '    [D]',
    '[N] [C]',
    '[Z] [M] [P]',
    ' 1   2   3',
    '',
    'move 1 from 2 to 1',
    'move 3 from 1 to 3',
    'move 2 from 2 to 1',
    'move 1 from 1 to 2',
].join('\r\n')

const processStackLine = (line) => {
    const output = []
    while (line) {
        const chunk = line.slice(0, 4)
        const data = chunk.replace(/[\[\]]/g, '').trim() || null
        output.push(data)
        line = line.slice(4)
    }
    return output
}


const processInput = (input, reverseOrder = true) => {
    let stack = input[0].split(/\r?\n/).map(line => processStackLine(line))
    // remove the footer
    stack.pop()
    if (reverseOrder) {
        // reverse the stack so the bottom is at the top
        stack = stack.reverse()
    }
    // less rearrange the data into columns
    const columns = []
    stack.forEach((stack, index) => {
        stack.forEach((val, i) => {
            columns[i] = columns[i] || []
            // if there is a value, add it to the column - otherwise ignore placeholders
            if (val) {
                columns[i].push(val)
            }
        })
    })

    const moves = input[1].trim().split(/\r?\n/).map(line => {
        const regex = /move (\d+) from (\d+) to (\d+)/
        const [, move, from, to] = line.match(regex)
        return {move: parseInt(move), from: parseInt(from), to: parseInt(to)}
    })
    return {columns, moves}
}

const executeMove = ({move, from, to, columns, reverse = true}) => {
    let moving = columns[from - 1].slice(-move)
    if (reverse) {
        moving = moving.reverse()
    }
    columns[from - 1] = columns[from - 1].slice(0, -move)
    columns[to - 1] = columns[to - 1].concat(moving)
}
const executeMoves = (moves, columns, reverse) => {
    moves.forEach(move => executeMove({...move, columns, reverse}))
}

const getTops = (columns) => columns.map(col => col[col.length - 1]).join('')

const main = async () => {
    const data = await fs.readFile(path.resolve(__dirname, 'input.txt'), 'utf8')
    const input = data.split(/\r?\n\r?\n/)
    const testData = testInput.split(/\r?\n\r?\n/)


    let {columns: testColumns, moves: testMoves} = processInput(testData)
    executeMoves(testMoves, testColumns)
    assert(getTops(testColumns) === 'CMZ', 'getTops(testColumns) failed')
    testColumns = processInput(testData).columns
    executeMoves(testMoves, testColumns, false)
    assert(getTops(testColumns) === 'MCD', `getTops(testColumnsCopy) failed - expected MCD got ${getTops(testColumns)}`)

    let {columns, moves} = processInput(input)
    executeMoves(moves, columns)
    console.log(`Part 1 Result: ${getTops(columns)}`)

    columns = processInput(input).columns
    executeMoves(moves, columns, false)
    console.log(`Part 2 Result: ${getTops(columns)}`)

    // console.log(`Part 2 Result: ${findPartialOverlaps(input)}`)


}

main().catch(console.error)
