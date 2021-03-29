if(block_count != 0){
    diff_bits = np.insert(diff_bits, 0, prev_sync_bits, axis=0);
}

position = 0;
while(True){
    block = diff_bits[position:position+26];
    potential_syndrome.resize(10,0);
    for(int i = 0 ; i<potential_syndrome.size() ; i++){
        for(int j = 0 ; j<26 ; j++){
            mult = block[j] && H[j,i];
            potential_syndrome[i] = (potential_syndrome[i] && !mult) || (!potential_syndrome[i] && mult);
        }
    }
    //convert to int
    potential_syndrome = static_cast<short int>(potential_syndrome);
    //Checks if syndrome A
    if((potential_syndrome).tolist() == [1,1,1,1,0,1,1,0,0,0]){
        if(last_position == -1 or printposition-last_position == 26){ 
            last_position = printposition;
            std::cerr << "Syndrome A at position " << printposition << std::endl;
            last_position = printposition;
        }
        else{
            std::cerr << "False positive Syndrome A at position " << printposition << std::endl;
        }
    }
    //Checks if syndrome B
    elif((potential_syndrome).tolist() == [1,1,1,1,0,1,0,1,0,0]){
        if(last_position == -1 or printposition-last_position == 26){ 
            std::cerr << "Syndrome B at position " << printposition << std::endl;
            last_position = printposition;
        }
        else{
            std::cerr << "False positive Syndrome B at position " << printposition << std::endl;
        }
    }
    //Checks if syndrome C
    elif((potential_syndrome).tolist() == [1,0,0,1,0,1,1,1,0,0]){
        if(last_position == -1 or printposition-last_position == 26){ 
            std::cerr << "Syndrome C at position " << printposition << std::endl;
            last_position = printposition;
        }
        else{
            std::cerr << "False positive Syndrome C at position " << printposition << std::endl;
        }
    }
    //Checks if syndrome D
    elif((potential_syndrome).tolist() == [1,0,0,1,0,1,1,0,0,0]){
        if(last_position == -1 or printposition-last_position == 26){ 
            std::cerr << "Syndrome D at position " << printposition << std::endl;
            last_position = printposition;
        }
        else{
            std::cerr << "False positive Syndrome D at position " << printposition << std::endl;
        }
    }
    //Breaks once it reaches the end
    position += 1;
    if(position+26 > diff_bits.size()-1){
        break;
    }
    printposition+=1;
}
//Creates list of bits not used 
prev_sync_bits = diff_bits[position-1::];

//Iterates through the blocks 
block_count += 1;


