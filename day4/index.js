const readInput = require('../libs/readInput')
const assert = require("assert")

const testInput = [
    '2-4,6-8',
    '2-3,4-5',
    '5-7,7-9',
    '2-8,3-7',
    '6-6,4-6',
    '2-6,4-8',
]

/**
 * A simple class to represent a pair of numbers
 */
class Pair {
    /**
     * @constructor
     * @param {string|number} x
     * @param {string|number} y
     */
    constructor(x, y) {
        this.x = parseInt(x)
        this.y = parseInt(y)
    }

    /**
     * Returns true if this pair fully overlaps with another pair
     * @param {Pair} other
     * @returns {boolean}
     */
    overlaps(other) {
        return (this.x <= other.x && this.y >= other.y) || (other.x <= this.x && other.y >= this.y)
    }

    /**
     * Returns true if this pair partially overlaps with another pair
     * @param {Pair} other
     * @returns {boolean}
     */
    partialOverlap(other) {
        return this.x <= other.y && this.y >= other.x
    }

    /**
     * Returns a new pair that is parsed from a string range (e.g. '1-2')
     * @param {string} range
     * @returns {Pair}
     */
    static fromRange(range) {
        const [a, b] = range.split('-').map(Number)
        return new Pair(a, b)
    }
}

const findOverlapsV2 = (input) => {
    const overlaps = {full: 0, partial: 0}
    input.forEach((line) => {
        const ranges = line.split(',') //.map((range) => Pair.fromRange(range))
        const a = Pair.fromRange(ranges[0])
        const b = Pair.fromRange(ranges[1])
        overlaps.full += ~~a.overlaps(b)
        overlaps.partial += ~~a.partialOverlap(b)
    })
    return overlaps
}

/**
 * Split lines into pairs then into ranges
 * @param {[string]} input
 * @returns {[[[number, number],[number, number]]]}
 */
const formatInput = (input) => input.map((line) => line.split(',').map((range) => range.split('-').map((n) => parseInt(n))))

// overlaps are when x1 <= x2 & y1 <= y2 or the reverse of pairs a and b
/**
 * Returns true if the two ranges overlap
 * overlaps are when x1 <= x2 & y1 <= y2 or the reverse of pairs a and b
 * @param {[number, number]} a
 * @param {[number, number]} b
 * @returns {boolean}
 */
const doesArrayOverlap = (a, b) => (a[0] <= b[0] && a[1] >= b[1]) || (a[0] >= b[0] && a[1] <= b[1])

const doesArrayPartiallyOverlap = (a, b) => {
    // do we have a full overlap?
    if (doesArrayOverlap(a, b)) return true
    // does a start on or before b ends and end on or after b starts
    if (a[0] <= b[0] && a[1] >= b[0]) return true
    // does b start on or before a starts and end on or before a ends
    return b[0] <= a[0] && b[1] >= a[0];
}

/**
 * Returns the number of overlaps. We do this by converting the true/false array into a number array using the ~ operator as in ~~true === 1 and ~~false === 0
 * @param {[boolean]} arr
 * @returns {number}
 */
const countOverlaps = (arr) => arr.reduce((acc, cur) => (acc += ~~cur, acc), 0)

/**
 * Takes an input array of strings and returns the number of overlaps
 * input -> formatInput() -> map() -> doesArrayOverlap() -> countOverlaps()
 * @param {[string]} input
 * @returns {number}
 */
const findOverlaps = (input) => countOverlaps(formatInput(input).map((arr) => doesArrayOverlap(...arr)))

/**
 * Takes an input array of strings and returns the number of partial overlaps
 * input -> formatInput() -> map() -> doesArrayPartiallyOverlap() -> countOverlaps()
 * @param {[string]} input
 * @returns {number}
 */
const findPartialOverlaps = (input) => countOverlaps(formatInput(input).map((arr) => doesArrayPartiallyOverlap(...arr)))

const main = async () => {
    const input = await readInput(__dirname, 'input.txt')

    assert(findOverlaps(testInput) === 2, 'findOverlaps failed')
    assert(findPartialOverlaps(testInput) === 4, 'findPartialOverlaps failed')
    assert(formatInput(testInput).map((arr) => doesArrayPartiallyOverlap(...arr)).join(',') === 'false,false,true,true,true,true', 'doesArrayPartiallyOverlap failed')

    console.log(`Part 1 Result: ${findOverlaps(input)}`)
    console.log(`Part 2 Result: ${findPartialOverlaps(input)}`)

    const testOverlaps = findOverlapsV2(testInput)
    const overlaps = findOverlapsV2(input)
    assert(testOverlaps.full === 2, 'findOverlapsV2.full failed')
    assert(testOverlaps.partial === 4, 'findOverlapsV2.partial failed')
    console.log(`More performant Part 1 Result: ${overlaps.full}`)
    console.log(`More performant Part 2 Result: ${overlaps.partial}`)
}

main().catch(console.error)
