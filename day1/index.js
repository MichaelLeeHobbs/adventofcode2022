const fs = require('fs/promises')
const path = require('path')

const main = async () => {
    const data = await fs.readFile(path.resolve(__dirname, 'input.txt'), 'utf8')
    const lines = data.split(/\r?\n/)
    const result = lines.reduce((acc, line) => {
        if (line === '') {
            acc.push(0)
            return acc
        }
        acc[acc.length - 1] += parseInt(line)
        return acc
    }, [0]).sort((a, b) => b - a)
    console.log(`Elf with most calories: ${result[0]}`)
    console.log(`Top 3 elves: ${result.slice(0, 3).reduce((acc, val) => acc + val, 0)}`)
}

main().catch(console.error)


const test = async () => {
    const data = await fs.readFile(path.resolve(__dirname, 'input.txt'), 'utf8')
    const result = data.split(/\r?\n\r?\n/)
        .map(line => line.split(/\r?\n/).reduce((acc, line) => acc + parseInt(line), 0))
        .sort((a, b) => b - a)
    console.log(`Elf with most calories: ${result[0]}`)
    console.log(`Top 3 elves: ${result.slice(0, 3).reduce((acc, val) => acc + val, 0)}`)
}
test().catch(console.error)
