import Foundation
import Metal

func readExactly(_ count: Int) throws -> Data? {
    var data=Data();data.reserveCapacity(count)
    while data.count<count {
        guard let block=try FileHandle.standardInput.read(upToCount:count-data.count), !block.isEmpty else {
            if data.isEmpty{return nil}
            throw NSError(domain:"SiderAtlas",code:1,userInfo:[NSLocalizedDescriptionKey:"Truncated input"])
        }
        data.append(block)
    }
    return data
}
func words(_ data: Data) -> [UInt32] {
    return data.withUnsafeBytes { ptr in
        (0..<data.count/4).map{UInt32(littleEndian:ptr.loadUnaligned(fromByteOffset:$0*4,as:UInt32.self))}
    }
}
do {
    guard let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue() else {
        throw NSError(domain:"SiderAtlas",code:2,userInfo:[NSLocalizedDescriptionKey:"GPU unavailable in this process"])
    }
    let source=try String(contentsOfFile:CommandLine.arguments[1],encoding:.utf8)
    let options=MTLCompileOptions();options.fastMathEnabled=false
    let library=try device.makeLibrary(source:source,options:options)
    guard let function=library.makeFunction(name:"project_material") else {fatalError("Missing kernel")}
    let pipeline=try device.makeComputePipelineState(function:function)
    let handshake=try JSONSerialization.data(withJSONObject:["device":device.name,"ready":true],options:[.sortedKeys])
    FileHandle.standardOutput.write(handshake);FileHandle.standardOutput.write(Data([10]))
    while let header=try readExactly(48) {
        try autoreleasepool {
            let h=words(header);guard h[0]==0x53414e36 else {fatalError("Wrong packet")}
            var buffers:[MTLBuffer]=[]
            for count in [h[8],h[10],h[9],h[11]] {
                guard let data=try readExactly(Int(count)) else {fatalError("Missing payload")}
                let buffer=data.withUnsafeBytes {p in device.makeBuffer(bytes:p.baseAddress!,length:data.count,options:.storageModeShared)!}
                buffers.append(buffer)
            }
            let dimensions=[h[1],h[2],h[3],h[4],h[5],h[6],h[7]]
            let dims=dimensions.withUnsafeBytes {device.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)!}
            let size=Int(h[1])*Int(h[2])*4*4
            let output=device.makeBuffer(length:size,options:.storageModeShared)!
            guard let command=queue.makeCommandBuffer(),let encoder=command.makeComputeCommandEncoder() else {fatalError("GPU encoder")}
            encoder.setComputePipelineState(pipeline)
            for i in 0..<4{encoder.setBuffer(buffers[i],offset:0,index:i)}
            encoder.setBuffer(dims,offset:0,index:4);encoder.setBuffer(output,offset:0,index:5)
            encoder.dispatchThreads(MTLSize(width:Int(h[1]),height:Int(h[2]),depth:1),
                                    threadsPerThreadgroup:MTLSize(width:16,height:8,depth:1))
            encoder.endEncoding();command.commit();command.waitUntilCompleted()
            if let error=command.error {throw error}
            let result:[UInt32]=[0x53414e36,UInt32(size)]
            result.withUnsafeBytes {FileHandle.standardOutput.write(Data($0))}
            FileHandle.standardOutput.write(Data(bytesNoCopy:output.contents(),count:size,deallocator:.none))
        }
    }
} catch {
    FileHandle.standardError.write(Data(("SiderAtlas GPU error: \(error)\n").utf8))
    exit(1)
}
