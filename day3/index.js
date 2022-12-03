const readInput = require('../libs/readInput')
const assert = require("assert")

const testInput = [
    'vJrwpWtwJgWrhcsFMMfFFhFp',
    'jqHRNqRjqzjGDLGLrsFMfFZSrLrFZsSL',
    'PmmdzqPrVvPwwTWBwg',
    'wMqvLMZHhHMvwLHjbvcjnnSBnvTQFn',
    'ttgJtRGJQctTZtZT',
    'CrZsJsPPZsGzwwsLwLmpwMDw',
]

const priority = (c) => {
    assert(c.length === 1, 'c must be a single character')
    const n = c.charCodeAt(0)
    return n > 96 ? n - 96 : n - 38
}


const splitter = (line) => {
    const splitIndex = line.length / 2
    return [line.slice(0, splitIndex).split(''), line.slice(splitIndex).split('')]
}

const getPriority = ([a, b]) => {
    let _priority = 0
    a.some((c) => {
        if (b.includes(c)) {
            _priority = priority(c)
            return true
        }
    })
    return _priority
}

const getBadge = ([a, b, c]) => {
    let badge = 0
    a.some((n) => {
        if (b.includes(n) && c.includes(n)) {
            badge = priority(n)
            return true
        }
    })
    return badge
}

/**
 * Accepts an array of strings and returns an array of arrays of strings with each child array containing three strings
 * @param {[string]} lines
 * @returns {[[string]]}
 */
const bucketer = (lines) => {
    return lines.reduce((acc, line) => {
        let bucketIndex = acc.length - 1
        const isFull = acc[bucketIndex].length === 3
        if (isFull) {
            acc.push([])
            bucketIndex++
        }
        acc[bucketIndex].push(line.split(''))
        return acc
    }, [[]])
}

const getPrioritySum = (input) => input.map(splitter).reduce((acc, [a, b]) => acc + getPriority([a, b]), 0)

const getBadgeSum = (input) => bucketer(input).reduce((acc, [a, b, c]) => acc + getBadge([a, b, c]), 0)

const main = async () => {
    const input = await readInput(__dirname, 'input.txt')

    assert(priority('a') === 1, 'a should have priority 1')
    assert(priority('A') === 27, 'A should have priority 27')
    assert(priority('z') === 26, 'z should have priority 26')
    assert(priority('Z') === 52, 'Z should have priority 52')
    assert(getPrioritySum(testInput) === 157, 'getPrioritySum failed')
    assert(getBadgeSum(testInput) === 70, 'getBadgeSum failed')

    console.log(`Part 1 Result: ${getPrioritySum(input)}`)
    console.log(`Part 2 Result: ${getBadgeSum(input)}`)
}

main().catch(console.error)
